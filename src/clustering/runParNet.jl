#addprocs(2)

addprocs(parse(Int,ARGS[6]))
#
input1 = ARGS[1]
input2 = ARGS[2]
OGGmap = ARGS[3]
@eval @everywhere output1 = $ARGS[4]
@eval @everywhere output2 = $ARGS[5]
boots = parse(Int,ARGS[7])
#
# #
# input2 = "../../../../Data/devnetwork/beeSmall"
# input1 = "../../../../Data/devnetwork/antSmall"
# OGGmap = "../../../../Data/devnetwork/OGGmap.txt"
# @everywhere output1 = "../../../../Data/devnetwork/antTest"
# @everywhere output2 = "../../../../Data/devnetwork/beeTest"
@everywhere include("smallNet_2species.jl")

#boots=3

pos1 = string(input1,"pos.txt")
neg1 = string(input1,"neg.txt")
pos2 = string(input2,"pos.txt")
neg2 = string(input2,"neg.txt")


@everywhere initial_temp = 1 #0.8 appears to be the right place to start when not considering orthologs
@everywhere epochs = 50000000 #taking ~4 hours on a single core with initial_temp = 0.8
@everywhere max_mods = 250
@everywhere coupling_constant = 3
@everywhere cooling_constant = 0.9

nGene1 = getNgene(pos1)
nGene1n = getNgene(neg1)
nGene1 = max(nGene1,nGene1n)
nGene2 = getNgene(pos2)
nGene2n = getNgene(neg2)
nGene2 = max(nGene2,nGene2n)
Adj1pos = getAdj(pos1,nGene1)
Adj1neg = getAdj(neg1,nGene1)
Adj2pos = getAdj(pos2,nGene2)
Adj2neg = getAdj(neg2,nGene2)
@eval @everywhere nGene = [$nGene1,$nGene2]


#Load orthology, make a dictionary
dict1,dict_ogg1,weights = genDictionary(OGGmap)
AdjO = OGGmat(dict1,dict_ogg1,weights)
@eval @everywhere dict = $dict1
@eval @everywhere dict_ogg = $dict_ogg1
@eval @everywhere AdjMatOGG = $AdjO

Adj_pos = [diagMat(M) for M in [Adj1pos,Adj2pos]]
Adj_neg = [diagMat(M) for M in [Adj1neg,Adj2neg]]

tot_pos = [sum(Adj_pos[j])/2 for j=1:2]
pos = [sum(Adj_pos[j],2) for j=1:2]
pos_adj = [(pos[i]*pos[i]')/(2*tot_pos[i]) for i=1:2]

tot_neg = [sum(Adj_neg[j])/2 for j=1:2]
neg = [sum(Adj_neg[j],2) for j=1:2]
neg_adj = [(neg[i]*neg[i]')/(2*tot_neg[i]) for i=1:2]

#We are basically just going to be evaluating the contribution of a pair of nodes to the network, so may as well
#compute the whole thing before hand
AdjMat2 = [Adj_pos[i]-Adj_neg[i]-pos_adj[i]+neg_adj[i] for i=1:2]
for n=1:2
    for i=1:nGene[n],j=1:nGene[n]
        if i == j
            AdjMat2[n][i,j] = 0
        end
    end
end

@eval @everywhere AdjMat = $AdjMat2

println("starting simulations for $boots times")

pmap((args)->runSim(args...),[rn for rn=1:boots])
