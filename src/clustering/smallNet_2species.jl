using Distributions
#using DistributedArrays

# addprocs(parse(Int,ARGS[4]))
addprocs(2)
#
# input1 = ARGS[1]
# input2 = ARGS[2]
# output1 = ARGS[3]
# output2 = ARGS[4]
input1 = "../../../../Data/devnetwork/beeNetSmall.csv"
input2 = "../../../../Data/devnetwork/ant100_net.csv"
output1 = "../../../../Data/devnetwork/beeNetSmall_clust.csv"
output2 = "../../../../Data/devnetwork/ant100_net_clust.csv"

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

temp = 0.1
epochs = 5
max_mods = 250
srand()

#Functions to update temperature and number of iterations per epoch
tf(t) = 0.7*t
itf(length) = floor(Int,1.0*length)

nGene1 = getNgene(input1)
nGene2 = getNgene(input2)
Adj1 = getAdj(input1,nGene1)
Adj2 = getAdj(input2,nGene2)

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
