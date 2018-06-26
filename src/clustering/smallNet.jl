using Distributions
#using DistributedArrays

inputfile = ARGS[1]
outputfile = ARGS[2]
nGene = parse(Int,ARGS[3])+1
println(nGene)
println("start")
spins = Vector{Int}(nGene)
test_spins = Vector{Int}(nGene)
Adj = falses(nGene,nGene)
temp = 0.1
epochs = 100
current_energy = Float64
MersenneTwister(0)
tf(t) = 0.5*t
itf(length) = floor(Int,1.0*length)

open(inputfile) do file
    for ln in eachline(file)
        sp = split(ln,"\t")
        a1 = parse(Int,sp[1]) + 1
        a2 = parse(Int,sp[2]) + 1
        Adj[a1,a2] = true
    end
end

tot_pos = sum(Adj)
pos_each = sum(Adj,1)

function CalcEnergy(spins)
    energy = 0
    for i=1:nGene,j=1:nGene
        if spins[i] == spins[j]
            energy = energy + Adj[i,j] - pos_each[i]*pos_each[j]/(2*tot_pos)
        end
    end
    #energy = -energy
    return energy
end

function Initialize()
    for i=1:nGene
        spins[i] = rand(1:nGene)
    end
    return spins
end

function move()
    node = rand(1:100)
    spins[node]=rand(1:nGene)
    return(spins)
end

function testMove(test_spins,true_spins)
    passed = 0
    testEnergy = CalcEnergy(test_spins)
    if testEnergy < current_energy
        true_spins = test_spins
        passed = 1
    else
        prob = exp((current_energy - testEnergy)/(temp))
        if prob > rand(Uniform(0,1),1)[1]
            true_spins = test_spins
            passed = 1
        end
    end
    return true_spins,passed
end

spins = Initialize()
initial = CalcEnergy(spins)

for iter=1:100
    success = 0
    for e=1:epochs
        current_energy = CalcEnergy(spins)
        test_spins = move()
        spins,passed = testMove(test_spins,spins)
        success = success+passed
    end
    println(success/epochs)
    temp = tf(temp)
    epochs = itf(epochs)
end

writedlm(outputfile,spins)
