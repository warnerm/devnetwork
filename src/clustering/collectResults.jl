addprocs(parse(Int,ARGS[4]))

#Identify modules as groups of genes that are in the same module 95% of the time
@eval @everywhere prefix1 = $ARGS[1]
@eval @everywhere prefix2 = $ARGS[2]
nFile = parse(Int,ARGS[3])

# prefix1 = "../../cluster_results/antOrthoNet"
# prefix2 = "../../cluster_results/beeOrthoNet"
# nFile = 3
function getNgene(file)
    lastGene = Int
    open(file) do f
        for ln in eachline(file)
            sp = split(ln,"\t")
            a1 = parse(Int,sp[1])
            a2 = parse(Int,sp[2])
            lastGene = a1
        end
    end
    return lastGene
end

nGene1 = getNgene(string(prefix1,"_1.txt"))
nGene2 = getNgene(string(prefix2,"_1.txt"))
@eval @everywhere nGene = [$nGene1, $nGene2]


@everywhere function AddConns(f)
    println(f)
    Aii = zeros(nGene[1],nGene[1])
    Aij = zeros(nGene[1],nGene[2])
    Ajj = zeros(nGene[2],nGene[2])
    f1 = open(readdlm,string(prefix1,"_",f,".txt"))
    f2 = open(readdlm,string(prefix2,"_",f,".txt"))
    for i=1:nGene[1],j=1:nGene[1]
        Aii[i,j] = (f1[i,2] == f1[j,2] ? Aii[i,j]+1 : Aii[i,j])
    end
    for i=1:nGene[2],j=1:nGene[2]
        Ajj[i,j] = (f2[i,2] == f2[j,2] ? Ajj[i,j]+1 : Ajj[i,j])
    end
    for i=1:nGene[1],j=1:nGene[2]
        Aij[i,j] = (f1[i,2] == f2[j,2] ? Aij[i,j]+1 : Aij[i,j])
    end
    return Array[Aii, Aij, Ajj]
end

AdjAll = pmap((args)->AddConns(args...),[f for f=1:nFile])
Adj = [zeros(nGene[1],nGene[1]),zeros(nGene[1],nGene[2]),zeros(nGene[2],nGene[2])]
for f=1:nFile
    for n=1:3
        Adj[n] = Adj[n] + AdjAll[f][n]
    end
end

Adj_sparse1 = [[1,1,i,j,Adj[1][i,j]/nFile] for i=1:nGene[1],j=1:nGene[2] if Adj[1][i,j] != 0 && i!=j]
Adj_sparse2 = [[1,2,i,j+nGene[1],Adj[2][i,j]/nFile] for i=1:nGene[1],j=1:nGene[2] if Adj[2][i,j] != 0 && i!=j]
Adj_sparse3 = [[2,2,i+nGene[1],j+nGene[1],Adj[3][i,j]/nFile] for i=1:nGene[2],j=1:nGene[2] if Adj[3][i,j] != 0 && i!=j]
Adj_sparse = [Adj_sparse1;Adj_sparse2;Adj_sparse3]

writedlm(string(prefix1,"_all_collectedRes.txt"),Adj_sparse)
