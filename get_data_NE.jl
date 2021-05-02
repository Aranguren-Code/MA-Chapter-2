# Data for the problem
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
using DataFrames
using CSV

Network = CSV.read("network_parameters.csv"); #change depending on depot setting
ArcsT1 = CSV.read("parameters_in_T1.csv");
ArcsT2 = CSV.read("parameters_in_T2.csv");
ArcsT3 = CSV.read("parameters_in_T3.csv");
NodesP = CSV.read("parameters_in_P.csv"); #rent included per hectare in column 6
NodesC = CSV.read("parameters_in_C.csv");
NodesD = CSV.read("parameters_in_D.csv"); 
Quality = CSV.read("P_Moist_Ash_Yield.csv")

# NodesD25 = CSV.read("parameters_in_D(2.5).csv");
# NodesD375 = CSV.read("parameters_in_D(3.75).csv");
# NodesD5 = CSV.read("parameters_in_D(5).csv")
Parameters = CSV.read("other_parameters.csv");
RCPs = CSV.read("RCP_Stochastic_Scenarios.csv");
Depots_SA = CSV.read("FINAL_SA_Solution.csv");

P = Network[1, 1] # total number of parcels
D = Network[1, 2] # total candidate depots
C = Network[1, 3] # total candidate plants
S = Network[1, 4] # total scenarios

NoArcsT1 = P * D
NoArcsT2 = D * C
NoArcsT3 = P * C
#Zeros arrays are standard Float64
TransCost1 = zeros(P, D)
TransCost2 = zeros(D, C)
TransCost3 = zeros(P, C)
D_Selected_Sol = zeros(D)
HarvCost = zeros(P)
Supply = zeros(P, S)
DepotCost = zeros(D)
DepotCap = zeros(D)
Demand = zeros(C)
ShortCost = zeros(C)
prob = zeros(S)
D_Model = zeros(Int8, D)

RentCost1 = Parameters[1, 2]
RentCost2 = Parameters[2, 2]
Bailing = Parameters[3, 2]
Fertilizer = Parameters[4, 2]
Swathing = Parameters[5, 2]
ParcelSize = Parameters[6, 2]
AvgKmph1 = Parameters[7, 2]
AvgKmph2 = Parameters[8, 2]
AvgKmph3 = Parameters[9, 2]
HrRate1 = Parameters[10, 2]
HrRate2 = Parameters[11, 2]
HrRate3 = Parameters[12, 2]
DryMass1 = Parameters[13, 2]
DryMass2 = Parameters[14, 2]
DryMass3 = Parameters[15, 2]
Load = Parameters[16, 2]
Unload = Parameters[17, 2]
Stack = Parameters[18, 2]

# transportation cost within the network via truck form parcels to depots (T1)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Demand From Coal Power plants
for i in 1:C 
Demand[i] = NodesC[i,3]
end

#Shortcost per power plant
for i in 1:C 
ShortCost[i] = NodesC[i,4]
end

i = 1 
j = 1
cnt = 1 

while cnt <= NoArcsT1
    if j <= D
        TransCost1[i, j] = 2 * ArcsT1[cnt, 4] / AvgKmph1 * HrRate1 / DryMass1
        TransCost1[i, j] = TransCost1[i, j] + Load + Unload + Stack
        global j += 1
        global cnt += 1
    else
        global i += 1
        global j = 1
    end
end

# transportation cost within the network via truck form depots to coal plants
#(T2)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

i = 1
j = 1
cnt = 1

while cnt <= NoArcsT2
    if j <= C
        TransCost2[i, j] = 2 * ArcsT2[cnt, 4] / AvgKmph2 * HrRate2 / DryMass2
        #TransCost2[i, j] = TransCost2[i, j] + Load + Unload + Stack
        #TransCost2[i, j] = TransCost2[i, j]  / 6
        #TransCost2[i, j] = TransCost2[i, j]  / (DryMass2 / DryMass1)
        global j += 1
        global cnt += 1
    else
        global i += 1
        global j = 1
    end
end

# transportation cost within the network via truck form parcels to coal plants
#(T3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

i = 1
j = 1
cnt = 1

while cnt <= NoArcsT3
    if j <= C
        TransCost3[i, j] = 2 * ArcsT3[cnt, 4] / AvgKmph3 * HrRate3 / DryMass3
        TransCost3[i, j] = TransCost3[i, j] + Load + Unload + Stack
        global j += 1
        global cnt += 1
    else
        global i += 1
        global j = 1
    end
end

# harvesting cost for every parcel in P
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for i in 1:P
    HarvCost[i] = Bailing + (Fertilizer + Swathing + NodesP[i,6]) / 
    (if NodesP[i,3] == 0.0 
        0.00000001
    else 
    NodesP[i, 3]
    end)
end

# supply of biomass in every parcel
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for i in 1:P
	for s in 1:S
		Supply[i,s] = 2500*RCPs[i,s+1]
	end
end

# investment cost to open a depot
# storage capacity of a depot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for i in 1:D
    DepotCost[i] = NodesD[i, 3]
    DepotCap[i] = NodesD[i, 4]
end

# demand at each power plant
# penalty cost for slack of supply at each power plant
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#for i in 1:C
#    Demand[i] = NodesC[i, 3]
#    ShortCost[i] = NodesC[i, 4]
#end

# probability for each scenarios
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for s in 1:S
	prob[s] = 1/S
end

Ash = zeros(P)
Ash_cost = 28.86
for a in 1:P
    Ash[a] = (Quality[a,6]/100)*Ash_cost
end

Moist = zeros(P)
for m in 1:P
    Moist[m] = Quality[m,5]
end
#for i in 1:D
#	D_Selected_Sol[i] = Depots_SA[i,2]
#end

#May be removed if not using simulated annealing
sim_par = CSV.read("SA_parameters.csv")
temp = Int(sim_par[1,2]) ; # maximum temperature (about 20-30) was max_temp
min_temp = Int(sim_par[2, 2]) # minimum temperature (set to zero)
epochs = Int(sim_par[3, 2]) # number of epochs (around the same as max_temp)
current_energy = (maximum(ShortCost) * sum(Demand)) #third party supplier cost 
Best = current_energy #initial solution
next_energy = 0 #initialize variable
SupplyMg = 0
cumm = 0 
DemandMg = Demand
max_dep = Int(ceil(sum(Demand)/minimum(DepotCap)))

base = zeros(Int64, D, 2) #ID and state
base_accept = zeros(Int64, D, 2)
previous_base = zeros(Int64, D, 2)
Swap_base = zeros(Int64, D,2)
Memory = zeros(Int64, max_dep, 3) #may delete
Swap_reject = zeros(Int64, D)

#This base actually uses the Depot ID not a random index
#because the number of available depots is smaller not the entire
#dataset

for i in 1:D
    base[i,1] = NodesD[i,2] 
    base[i,2] = 0
end

for i in 1:D
    base_accept[i,1] = NodesD[i,2] 
    base_accept[i,2] = 0
end

# previous_base = base
temp = 10
swap_variable = 0
sum_base = 0
s1 = 0
s0 = 0 
b1 = 0
b0 = 0 
ite = 0
id1 = 0
op_id = 0
no_depots = 0
delta_d = 0
delta_e = 0
count = 0 
r_accept = 0 
flag = 0
Best_base = base

