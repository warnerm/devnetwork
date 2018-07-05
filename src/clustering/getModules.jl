# inputfile = "antOrthoNet_all_collectedRes.txt"
# outputfile = "test"
inputfile = ARGS[1]
outputfile = ARGS[2]

f = open(readdlm,inputfile)
genes_remaining = collect(1:Int(maximum(f[:,3:4]))) #make list of all genes

f = f[f[:,5] >= 0.8,:]
mod = 0
ModDef = Vector{Int}(length(genes_remaining))

while length(genes_remaining) > 0
    println(length(genes_remaining))
    mod = mod + 1
    new_gene = genes_remaining[1]
    splice!(genes_remaining,1)
    next_genes = []
    found = 0
    while true
        ModDef[new_gene] = mod
        conns = [Int(f[row,4]) for row in 1:size(f,1) if f[row,3] == new_gene]
        conns2 = [Int(f[row,3]) for row in 1:size(f,1) if f[row,4] == new_gene]
        append!(conns,conns2)
        if length(conns) > 0
            append!(next_genes,conns)
        end
        while true
            if length(next_genes)==0
                found = 0
                break
            end
            next = next_genes[1]
            if next in genes_remaining
                filter!(x->x!=next,genes_remaining)
                new_gene = next
                splice!(next_genes,1)
                found = 1
                break
            else
                splice!(next_genes,1)
            end
        end
        if found == 0
            break
        end
    end
end

sOut = [1:length(ModDef) ModDef]

writedlm(outputfile,sOut)
