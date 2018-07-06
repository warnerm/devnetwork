@everywhere using Distributions
#using DistributedArrays
#
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
begin
    @everywhere function Initialize(nGene)
        spin = Vector{Int}(nGene)
        for i=1:nGene
            spin[i] = rand(1:max_mods)
        end
        return spin
    end
end

#calculate the energy for a given node/spin, given the rest of the network as is
@everywhere function calcPartial(spins,node,indiv_spin,net)
    sameMod = [j for j in  1:nGene[net] if spins[net][j] == indiv_spin]
    if length(sameMod) > 0
        energy_all = sum([AdjMat[net][node,j] for j in sameMod])
    else
        energy_all = 0
    end
    if haskey(dict[net],node)
        OGG = dict[net][node]
        partners = dict_ogg[3-net][OGG]
        partners = [p for p in partners if p <= nGene[3-net]]
        sameMod_OGG = [j for j in partners if spins[3-net][j]==indiv_spin]
        if length(sameMod_OGG) > 0
            if net == 1
                energy_all = energy_all+sum([AdjMatOGG[node,j] for j in sameMod_OGG])
            else
                energy_all = energy_all+sum([AdjMatOGG[j,node] for j in sameMod_OGG])
            end
        end
    end
    return -energy_all
end

begin
    @everywhere function move(spins,temp)
        net = rand(1:2)
        node = rand(1:nGene[net])
        edit_spin = copy(spins[net][node])
        while edit_spin == spins[net][node]
            edit_spin = rand(1:max_mods)
        end
        oldEnergy = calcPartial(spins,node,spins[net][node],net)
        newEnergy = calcPartial(spins,node,edit_spin,net)
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

function OGGmat(dict,dict_ogg,weights)
    matO = zeros(nGene[1],nGene[2])
    for i=1:nGene[1]
        if haskey(dict[1],i)
            OGG = dict[1][i]
            partners = dict_ogg[2][OGG]
            if length(partners) > 0
                for j in partners
                    if j <= nGene[2]
                        matO[i,j] = weights[OGG]
                    end
                end
            end
        end
    end
    return matO
end

@everywhere function runSim(run)
    srand()
    spins = [Initialize(nGene[i]) for i=1:2]
    temp = initial_temp
    while true
        success = 0
        for e=1:epochs
            passed = move(spins,temp)
            success = success+passed
        end
        println(success/epochs)
        if success/epochs < 0.01 #Stop iterations if there are very few successes
            break
        end
        temp = temp*0.9
    end

    sOut1 = [1:nGene[1] spins[1]]
    sOut2 = [1:nGene[2] spins[2]]

    writedlm(string(output1,"Net_",run,".txt"),sOut1)
    writedlm(string(output2,"Net_",run,".txt"),sOut2)
    return 0
end
