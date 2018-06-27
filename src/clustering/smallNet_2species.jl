using Distributions
#using DistributedArrays

# addprocs(parse(Int,ARGS[4]))
addprocs(2)
#
# input1 = ARGS[1]
# input2 = ARGS[2]
# output1 = ARGS[3]
# output2 = ARGS[4]
input2 = "../../../../Data/devnetwork/beeNetSmall.csv"
input1 = "../../../../Data/devnetwork/ant100_net.csv"
output2 = "../../../../Data/devnetwork/beeNetSmall_clust.csv"
output1 = "../../../../Data/devnetwork/ant100_net_clust.csv"
OGGfile = "../../../../Data/devnetwork/OGGmap.txt"

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
    spin = Vector{Int64}(nGene)
    for i=1:nGene
        spin[i] = rand(1:max_mods)
    end
    return spin
end

#calculate the energy for a given node/spin, given the rest of the network as is
function calcPartial(node,indiv_spin,net)
    @eval @everywhere i = $node
    @eval @everywhere s = $indiv_spin
    @eval @everywhere n = $net
    energy = pmap((args)->nodeEnergy(args...),[j for j=1:nGene[net]])
    #Add energy if orthologs are in the same module
    ortholog = get(dict[net],node,-1)
    if ortholog != - 1 && ortholog <= nGene[3-net]
        if spins[3-net][ortholog] == spins[net][node]
            energy = energy + coupling_constant*weights[net][node]
        end
    end
    return -sum(energy)
end

@everywhere begin
    #Objective function to test energy with respect to given node
    function nodeEnergy(j)
        if spins[n][j] ==  s
            return ((Adj_all[n][i,j]+Adj_all[n][j,i]) -  2*pos_each[n][i]*pos_each[n][j]/(2*tot_pos[n]))
        else
            return 0
        end
    end
end

function move()
    net = rand(1:2)
    node = rand(1:nGene[net])
    edit_spin = rand(1:max_mods)
    oldEnergy = calcPartial(node,spins[net][node],net)
    newEnergy = calcPartial(node,edit_spin,net)
    passed = 0
    @eval @everywhere test_node = $node
    if newEnergy < oldEnergy
        @eval @everywhere spins[$net][test_node] = $edit_spin
        passed = 1
    else
        prob = exp((oldEnergy - newEnergy)/(temp))
        if prob > rand(Uniform(0,1),1)[1]
            @eval @everywhere spins[$net][test_node] = $edit_spin
            passed = 1
        end
    end
    return passed
end

function genDictionary(OGGfile)
    dict1 = Dict{Int64,Int64}()
    dict2 = Dict{Int64,Int64}()
    weights1 = Dict{Int64,Float64}()
    weights2 = Dict{Int64,Float64}()
    open(OGGfile) do f
        for ln in eachline(f)
            sp = split(ln,"\t")
            a1 = parse(Int,sp[1])
            a2 = parse(Int,sp[2])
            a3 = parse(Float,sp[3])
            dict1[a1] = a2
            dict2[a2] = a1
            weights1[a1] = a3
            weights2[a2] = a3

        end
    end
    all_dict = [dict1,dict2]
    weights = [weights1,weights2]
    return all_dict,weights
end

temp = 0.1
epochs = 5
max_mods = 250
coupling_constant = 3
srand()

#Functions to update temperature and number of iterations per epoch
tf(t) = 0.7*t
itf(length) = floor(Int,1.0*length)

nGene1 = getNgene(input1)
nGene2 = getNgene(input2)
Adj1 = getAdj(input1,nGene1)
Adj2 = getAdj(input2,nGene2)

#Load orthology, make a dictionary
dict,weights = genDictionary(OGGfile)


@eval @everywhere nGene = [$nGene1,$nGene2]
@eval @everywhere Adj_all = [$Adj1,$Adj2]

@everywhere tot_pos = [sum(Adj_all[j]) for j=1:2]
@everywhere pos_each = [sum(Adj_all[j],1) for j=1:2]
spin_list = [Initialize(nGene[i]) for i=1:2]
@eval @everywhere spins = $spin_list


for iter=1:2
    success = 0
    for e=1:epochs
        passed = move()
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
