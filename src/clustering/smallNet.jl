using Distributions
#using DistributedArrays

addprocs(parse(Int,ARGS[4]))

inputfile = ARGS[1]
outputfile = ARGS[2]
nGene = parse(Int,ARGS[3])+1
@eval @everywhere num = $nGene
println(nGene)
println("start")
test_spins = Vector{Int}(nGene)
Adj = falses(nGene,nGene)

@everywhere temp = 0.1
@everywhere epochs = 100
@everywhere current_energy = Float64
max_mods = 250
srand()
tf(t) = 0.7*t
itf(length) = floor(Int,1.0*length)

open(inputfile) do file
    for ln in eachline(file)
        sp = split(ln,"\t")
        a1 = parse(Int,sp[1]) + 1
        a2 = parse(Int,sp[2]) + 1
        Adj[a1,a2] = true
    end
end

@eval @everywhere Adj_all = $Adj

@everywhere tot_pos = sum(Adj_all)
@everywhere pos_each = sum(Adj_all,1)

@everywhere begin
    function TotEnergy(i,j)
        if spins[i] == spins[j]
            return Adj_all[i,j] -  pos_each[i]*pos_each[j]/(2*tot_pos)
        else
            return 0
        end
    end
end

function CalcEnergy()
    energy = pmap((args)->TotEnergy(args...),[[i,j] for i=1:num,j=1:num])
    return -sum(energy)
end

function Initialize()
    for i=1:nGene
        test_spins[i] = rand(1:max_mods)
    end
    return test_spins
end

@everywhere begin
    function nodeEnergy(j)
        if spins[j] ==  s
            return ((Adj_all[i,j]+Adj_all[j,i]) -  2*pos_each[i]*pos_each[j]/(2*tot_pos))
        else
            return 0
        end
    end
end

function calcPartial(node,indiv_spin)
    @eval @everywhere i = $node
    @eval @everywhere s = $indiv_spin
    energy = pmap((args)->nodeEnergy(args...),[j for j=1:num])
    return -sum(energy)
end

function move()
    node = rand(1:100)
    edit_spin = rand(1:max_mods)
    oldEnergy = calcPartial(node,spins[node])
    newEnergy = calcPartial(node,edit_spin)
    passed = 0
    @eval @everywhere test_node = $node
    if newEnergy < oldEnergy
        @eval @everywhere spins[test_node] = $edit_spin
        passed = 1
    else
        prob = exp((oldEnergy - newEnergy)/(temp))
        if prob > rand(Uniform(0,1),1)[1]
            @eval @everywhere spins[test_node] = $edit_spin
            passed = 1
        end
    end
    return passed
end

total_spins = Initialize()
@eval @everywhere spins = $total_spins
#
for iter=1:10000
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

writedlm(outputfile,spins)
