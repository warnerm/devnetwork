using Distributions
#using DistributedArrays
#
# addprocs(parse(Int,ARGS[4]))
#
# input1 = ARGS[1]
# input2 = ARGS[2]
# OGGmap = ARGS[3]


input2 = "../../../../Data/devnetwork/beeSmall"
input1 = "../../../../Data/devnetwork/antSmall"
OGGmap = "../../../../Data/devnetwork/OGGmap.txt"

pos1 = string(input1,"pos.txt")
neg1 = string(input1,"neg.txt")
pos2 = string(input2,"pos.txt")
neg2 = string(input2,"neg.txt")
output1 = string(input1,"net_clust.csv")
output2 = string(input2,"net_clust.csv")

#First, find number of genes
function getNgene(file)
    lastGene = Int
    open(file) do f
        for ln in eachline(file)
            sp = split(ln,"\t")
            a1 = parse(Int,sp[1]) + 1
            a2 = parse(Int,sp[2]) + 1
            lastGene = a1
        end
    end
    return lastGene
end

#From file giving list of directed connections, make an adjacency matrix
function getAdj(file,nGene)
    Adj = falses(nGene,nGene)
    open(file) do f
        for ln in eachline(file)
            sp = split(ln,"\t")
            a1 = parse(Int,sp[1]) + 1
            a2 = parse(Int,sp[2]) + 1
            Adj[a1,a2] = true
        end
    end
    return Adj
end

#Intialize spins with random modules
function Initialize(nGene)
    spin = Vector{UInt8}(nGene)
    for i=1:nGene
        spin[i] = convert(UInt8,rand(1:max_mods))
    end
    return spin
end

#calculate the energy for a given node/spin, given the rest of the network as is
function calcPartial(node,indiv_spin,net)
    @eval @everywhere i = $node
    @eval @everywhere n = $net
    sameMod = [j for j in  1:nGene[net] if spins[net][j] == indiv_spin]
    if length(sameMod) > 0
        energy_all = sum([AdjMat[net][node,j] for j in sameMod])
    else
        energy_all = 0
    end
    #Add energy if orthologs are in the same module
    if haskey(dict[net],node)
        OGG = dict[net][node]
        partners = dict_ogg[3-net][OGG]
        numSame = 0
        for i=1:length(partners)
            if partners[i] != 0 && partners[i] <= nGene[3-net]
                if spins[3-net][partners[i]] == indiv_spin
                    numSame = numSame + 1
                end
            end
        end
        energy_all = energy_all + coupling_constant*weights[OGG]*numSame
    end
    return -energy_all
end
#
# @everywhere begin
#     #Objective function to test energy with respect to given node
#     function nodeEnergy(j)
#         return (Adj_all_pos[n][i,j]-Adj_all_neg[n][i,j] - posA[n][i,j] + negA[n][i,j])
#     end
# end

function move()
    net = rand(1:2)
    node = rand(1:nGene[net])
    edit_spin = copy(spins[net][node])
    while edit_spin == spins[net][node]
        edit_spin = convert(UInt8,rand(1:max_mods))
    end
    oldEnergy = calcPartial(node,spins[net][node],net)
    newEnergy = calcPartial(node,edit_spin,net)
    passed = 0
    if newEnergy < oldEnergy
        spins[net][node] = edit_spin
        passed = 1
    else
        prob = exp(-(newEnergy - oldEnergy)/(temp))
        test_num = rand(Uniform(0,1))
        if prob > test_num
            spins[net][node] = edit_spin
            passed = 1
        end
    end
    return passed
end

function genDictionary(OGGfile)
    dict1 = Dict{Int64,String}()
    dict2 = Dict{Int64,String}()
    ogg1 = Dict{String,Vector{Int64}}()
    ogg2 = Dict{String,Vector{Int64}}()
    weights = Dict{String,Float64}()
    open(OGGfile) do f
        for ln in eachline(f)
            sp = split(ln,"\t")
            a1 = parse(Int,sp[1])+1
            a2 = parse(Int,sp[2])+1
            a3 = sp[3]
            a4 = parse(Float64,sp[4])
            dict1[a1] = a3
            dict2[a2] = a3
            if haskey(ogg1,a3)
                ogg1[a3] = push!(ogg1[a3],a1)
                ogg2[a3] = push!(ogg1[a3],a2)
            else
                ogg1[a3] = [a1]
                ogg2[a3] = [a2]
            end
            weights[a3] = a4
        end
    end
    all_dict = [dict1,dict2]
    all_ogg = [ogg1,ogg2]
    return all_dict,all_ogg,weights
end

function diagMat(M)
    newMat = [M[i,j] + M[j,i] for i=1:size(M,1),j=1:size(M,1)]
    return newMat
end

temp = .05
epochs = 10
max_mods = 250
coupling_constant = 3
srand()

#Functions to update temperature and number of iterations per epoch
tf(t) = 0.9*t
itf(length) = floor(Int,1.0*length)

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

#Load orthology, make a dictionary
dict,dict_ogg,weights = genDictionary(OGGmap)

nGene = [nGene1,nGene2]
Adj_pos = [diagMat(M) for M in [Adj1pos,Adj2pos]]
Adj_neg = [diagMat(M) for M in [Adj1neg,Adj2neg]]

tot_pos = [sum(Adj_all_pos[j]) for j=1:2]
pos = [sum(Adj_all_pos[j],2) for j=1:2]
pos_adj = [(pos[i]*pos[i]')/(2*tot_pos[i]) for i=1:2]

tot_neg = [sum(Adj_all_neg[j]) for j=1:2]
neg = [sum(Adj_all_neg[j],2) for j=1:2]
neg_adj = [(neg[i]*neg[i]')/(2*tot_neg[i]) for i=1:2]

#We are basically just going to be evaluating the contribution of a pair of nodes to the network, so may as well
#compute the whole thing before hand
AdjMat = [Adj_pos[i]-Adj_neg[i]-pos_adj[i]+neg_adj[i] for i=1:2]

spins = [Initialize(nGene[i]) for i=1:2]

for iter=1:10000
    success = 0
    for e=1:epochs
        @time passed = move()
        success = success+passed
    end
    println(success/epochs)
    if success/epochs < 0.01 #Stop iterations if there are very few successes
        break
    end
    temp = tf(temp)
    epochs = itf(epochs)
end

writedlm(output1,spins[1])
writedlm(output2,spins[2])
