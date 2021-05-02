#This program is version 2.14 of the simulated annealing

#The improvements involve making it much faster and 
# reduce computational burden 

# Attempt in this version: Remove Non-base from the code 

#Load packages

using JuMP
using CPLEX
using CSV
using DataFrames
#Load files 


#include("get_data_NE.jl")
include("Hub_Spoke_P3.jl")
sim_par = CSV.read("SA_parameters.csv");
temp = Int(sim_par[1,2]) ; # maximum temperature (about 20-30) was max_temp
min_temp = Int(sim_par[2, 2]) # minimum temperature (set to zero)
epochs = Int(sim_par[3, 2]) # number of epochs (around the same as max_temp)

# Initial Conditions 

current_energy = (maximum(ShortCost) * sum(Demand)) #third party supplier cost 
Best = current_energy #initial solution
next_energy = 0 #initialize variable
SupplyMg = 0
cumm = 0 
DemandMg = Demand
max_dep = Int(ceil(sum(Demand)/minimum(DepotCap)))

# Initialize data structures -----------------------------------------------------

base = zeros(Int8, D, 2)
Swap_base = zeros(Int8, D,2)
Memory = zeros(Int8, max_dep, 3) #may delete
Swap_reject = zeros(Int8, D)

for i in 1:D
    base[i,1] = (i)
    base[i,2] = (0)
end

# variables to use within while/for loop
sum_base = 0
ite = 0
op_id = 0
no_depots = 0
delta_d = 0
delta_e = 0
count = 0 
r_accept = 0 
flag = 0
Best_base = base

#Delete some of these later
best_output = open("Solution.csv", "w")
swaps = open("swaps.csv", "w")
all_base = open("All Bases.csv", "w")
f1 = open("checking.csv","w")
write(f1,"temp, ite, no_depots, delta_d, depot_id, status, op_id, \r\n")
write(best_output,"temp, epochs, move, id, depots_selected, Solution, \r\n")
write(swaps, "i, temp, epochs, id1, id2, \r\n")
write(all_base,"temp, epochs, id, depots_selected, delta_e, flag,\r\n")
start = time()  #know how long, put toc() at end

while temp >= min_temp && count <= 10

    ite = 0 
    while ite <= epochs 
        no_depots = Int(ceil(rand()*max_dep))
        
        delta_d = no_depots - sum_base #main determinant for move
        
        if delta_d > 0 #Addition
            op_id = 1

            for i in 1:delta_d
            r1 = Int(ceil(rand()*(D-sum_base)))
            id1 = base[r1,1]
            JuMP.fix(D_Model[id1], 1)
            Memory[i, 1] = r1
            Memory[i, 2] = id1
            Memory[i, 3] = op_id 
            end
        
        elseif delta_d < 0 #Removal
            op_id = 0

            delta_d = abs(delta_d)
            for i in 1:delta_d
            r1 = Int(ceil(rand()*sum_base)) + (D-sum_base)
            id1 = base[r1,1]
            JuMP.fix(D_Model[id1],0)
            Memory[i, 1] = r1
            Memory[i, 2] = id1
            Memory[i, 3] = op_id
            end
        
        else #SWAPPING
        op_id = 2
        Swap_base = base
            for i in 1:Int(ceil(D*exp(-temp/D*(factorial(max_dep)))))
                s1 = Int(ceil(rand()*D))
                s2 = Int(ceil(rand()*D))
                id1 = Swap_base[s1,1]
                id2 = Swap_base[s2,1]
                b1 = Swap_base[s1,2]
                b2 = Swap_base[s2,2]
                Swap_base[s1,2] = b2
                Swap_base[s2,2] = b1
                JuMP.fix(D_Model[id1],b2)
                JuMP.fix(D_Model[id2],b1)
                #CSV.write(swaps, "$i, $temp, $ite, $id1, $id2, \r\n") #delete later
            end
        end        

        if op_id == 2
            for i in 1:D
            write(f1,"$temp, $ite, $no_depots, $delta_d,")
            write(f1,"$(Swap_base[i,1]), $(Swap_base[i,2]), $op_id,\r\n")
            end
        else
            for i in 1:delta_d
            write(f1,"$temp, $ite, $no_depots, $delta_d,")
            write(f1,"$(Memory[i,2]), r1 = $(Memory[i,1]), $(Memory[i,3]),\r\n")
            end
        end

        #Obtained objective function value by solving prob:
        optimize!(HS_P3)
        global next_energy = JuMP.objective_value(HS_P3)
        global delta_e = current_energy - next_energy 

        #Delete Later - Flag for Neighborhood move
		if delta_e > 0
			flag = 1
		elseif delta_e < 0 
			if exp(delta_e/temp) > r_accept
				flag = 1
			else 
				flag = 0
			end
		else
			flag = 0
		end
        
        for i in 1:D
        write(all_base, "$temp,$ite,$(base[i,1]), $(base[i,2]), $delta_e, $flag, \r\n")
        end

        if delta_e > 0 #Accept all moves
        
            if op_id == 0 #Remove
                for i in 1:delta_d
                    r0 = Memory[i,1] #random number 
                    base[r0,2] = 0 #remove from base
                end
            elseif op_id == 1 #Addition
                for i in 1:delta_d
                    r1 = Memory[i,1]
                    base[r1,2] = 1 #add to base          
                end
            elseif op_id == 2 #SWAPPING
                base = Swap_base #accept the swap base
            end

            global current_energy = next_energy 
            global base = sortslices(base, dims = 1 , by = x -> x[2])
            global sum_base = sum(base[:,2])
            
            #Record Best Base and Reset count
            if next_energy < Best 
                global Best = next_energy 
                global Best_base = base
                count = 0 # Reset timer
                for i in 1:D
                    write(best_output, "$temp, $ite, $op_id, $(Best_base[i,1]),")
                    write(best_output, "$(Best_base[i,2]), $Best,\r\n" )
                end
            end 
        
        elseif delta_e < 0
            r_acccept = rand()
            if exp(delta_e/temp) > r_accept #if met we accept move
                
                if op_id == 0 
                    for i in 1:delta_d
                        id = Memory[i,1]
                        base[id,2] = 0 #remove from base
                    end                
                elseif op_id == 1 
                    for i in  1:delta_d
                        id = Memory[i,1]
                        base[id,2] = 1 #add to base          
                    end
                elseif op_id == 2 
                    base = Swap_base #accept the swap base
                end
            
            else  #metropolis criteria not met, movement denied

                if op_id == 0 #Removal
                    for i in 1:delta_d
                        id0 = Memory[i,2]
                        JuMP.fix(D_Model[id0], 1) #Add
                    end
                elseif op_id == 1 
                    for i in 1:delta_d
                        id1 = Memory[i,2]
                        JuMP.fix(D_Model[id1],0)
                    end                
                elseif op_id == 2
                    Swap_reject = Swap_base[:,2] - base[:,2]
                    for i in 1:D
                        if  Swap_reject[i,1] == -1 #means removed depot
                            id_swap = base[i,1]
                            JuMP.fix(D_Model[id_swap],1) #add it back
                        elseif Swap_reject[i,1] == 1 #means added depots
                            id_swap = base[i,1] 
                            JuMP.fix(D_Model[id_swap],0) #remove it again
                        end
                    end
                    Swap_base = base #reset swap_base
                end 

            end
            
            current_energy = next_energy
            global base = sortslices(base, dims = 1 , by = x -> x[2])
            sum_base = sum(base[:,2])

            #Record Best Base and reset count
            if next_energy < Best 
                Best = next_energy 
                Best_base = base
                count = 0 # Reset timer
                for i in 1:D
                    write(best_output, "$temp, $ite, $op_id, $(Best_base[i,1]),")
                    write(best_output, "$(Best_base[i,2]), $Best,\r\n" )
                end
            end 
        
        else # if the delta_e is zero it's not accepted

            if op_id == 0 #removed rejected
                for i in 1:delta_d
                    id0 = Memory[i,2]
                    JuMP.fix(D_Model[id0], 1) #Add back
                end
            elseif op_id == 1 
                for i in 1:delta_d
                    id1 = Memory[i,2]
                    JuMP.fix(D_Model[id1],0)
                end 
            elseif op_id == 2
                Swap_reject = Swap_base[:,2] - base[:,2]
                    for i in 1:D
                        if  Swap_reject[i,1] == -1 #means removed depot
                            id_swap = base[i,1]
                            JuMP.fix(D_Model[id_swap],1) #add it back
                        elseif Swap_reject[i,1] == 1 #means added depots
                            id_swap = base[i,1] 
                            JuMP.fix(D_Model[id_swap],0) #remove it again
                        end
                    end
                Swap_base = base #reset swap_base
            end
            
            #This is after the Base has been changed
            # if delta_d == 0 # SWAPPING
			#     for i in 1:(sum(Swap_base[:,2]))
			#     write(f1, "$temp,$ite,$no_depots,$delta_d,")
			#     write(f1, "$op_id,$(Swap_base[i, 1]),$(Swap_base[i, 2]),")
			#     write(f1, "$current_energy,$next_energy,$delta_e,$Best,$flag,\r\n")
			#     end
		    # else #for both deleting and adding
			#     for i in 1:delta_d
			#     write(f1, "$temp,$ite,$no_depots,$delta_d,")
			#     write(f1, "$op_id,$(Memory[i, 1]),$(Memory[i, 2]),")
			#     write(f1, "$current_energy,$next_energy,$delta_e,$Best,$flag,\r\n")
			# end
        end
    ite += 1
    end
    global temp -= 1  
    global count += 1 
end

close(swaps)
close(best_output)
close(f1)
#close(all_base)
elapse = time() - start