using Distributions
#using DistributedArrays

# addprocs(parse(Int,ARGS[6]))
# #
# input1 = ARGS[1]
# input2 = ARGS[2]
# output1 = ARGS[3]
# output2 = ARGS[4]
# OGGmap = ARGS[5]
input2 = "../../../../Data/devnetwork/beeNetSmall.csv"
input1 = "../../../../Data/devnetwork/ant100_net.csv"
output2 = "../../../../Data/devnetwork/beeNetSmall_clust.csv"
output1 = "../../../../Data/devnetwork/ant100_net_clust.csv"
OGGmap = "../../../../Data/devnetwork/OGGmap.txt"
addprocs(2)

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
    if length(sameMod) != 0
        energy = pmap((args)->nodeEnergy(args...),sameMod)
        energy_all = sum(energy)
    else
        energy_all = 0
    end
    #Add energy if orthologs are in the same module
    if haskey(dict[net],node)
        OGG = dict[net][node]
        partners = dict_ogg[3-net][OGG]
        for i=1:length(partners)
            if partners[i] != 0 && partners[i] <= nGene[3-net]
                if spins[3-net][partners[i]] == indiv_spin
                    energy_all = energy_all + coupling_constant*weights[OGG]
                end
            end
        end
    end
    return -energy_all
end

@everywhere begin
    #Objective function to test energy with respect to given node
    function nodeEnergy(j)
        return (Adj_all[n][i,j]+Adj_all[n][j,i] -  pos_out[n][i]*pos_in[n][j]/(2*tot_pos[n]) - pos_in[n][i]*pos_out[n][j]/(2*tot_pos[n]))
    end
end

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
    @eval @everywhere test_node = $node
    if newEnergy < oldEnergy
        spins[net][test_node] = edit_spin
        passed = 1
    else
        prob = exp(-(newEnergy - oldEnergy)/(temp))
        test_num = rand(Uniform(0,1))
        if prob > test_num
            spins[net][test_node] = edit_spin
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

temp = .05
epochs = 50
max_mods = 250
coupling_constant = 3
srand()

#Functions to update temperature and number of iterations per epoch
tf(t) = 0.5*t
itf(length) = floor(Int,1.0*length)

nGene1 = getNgene(input1)
nGene2 = getNgene(input2)
Adj1 = getAdj(input1,nGene1)
Adj2 = getAdj(input2,nGene2)

#Load orthology, make a dictionary
dict,dict_ogg,weights = genDictionary(OGGmap)


nGene = [nGene1,nGene2]
@eval @everywhere Adj_all = [$Adj1,$Adj2]

@everywhere tot_pos = [sum(Adj_all[j]) for j=1:2]
@everywhere pos_in = [sum(Adj_all[j],1) for j=1:2]
@everywhere pos_out = [sum(Adj_all[j],2) for j=1:2]
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
