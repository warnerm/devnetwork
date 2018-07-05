# inputfile = "../../../../Data/devnetwork/ant100_net.csv"
# adjFile = "../../../../Data/devnetwork/antsSMALLpCor.csv"
# outPre = "../../../../Data/devnetwork/antSmall"
#
inputfile = ARGS[1]
adjFile = ARGS[2]
outPre = ARGS[3]

outPos = string(outPre,"pos.txt")
outNeg = string(outPre,"neg.txt")
open(outNeg,"w") do out
end
open(outPos,"w") do out
end

open(inputfile) do f
    for line in eachline(f)
        sp = split(line,"\t")
        a1 = parse(Int,sp[1])+1
        a2 = parse(Int,sp[2])+1
        corList = open(readlines,adjFile)[a1+1]
        adj = split(corList,",")
        adjT = parse(Float64,adj[a2])
        if adjT < 0
            open(outNeg,"a") do out
                write(out,line,"\n")
            end
        else
            open(outPos,"a") do out
                write(out,line,"\n")
            end
        end
    end
end
