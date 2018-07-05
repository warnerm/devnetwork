inputfile = ARGS[1]
outputfile = ARGS[2]

f = open(readdlm,inputfile)
genes_remaining = [1:maximum(f[:,3:4])] #make list of all genes
mod = 0
ModDef = Vector{Int}(length(genes_remaining))

while length(genes_remaining) > 0
    println(length(genes_remaining))
    mod = mod + 1
    new_gene = genes_remaining[1]
    ModDef[new_gene] = mod
    next_genes = []
    found = 0
    while true
        conns = filter(x -> x[3] == new_gene && x[5] >= 0.95,f)
        if size(conns,1) > 0
            append!(next_genes,conns[:,4])
            [ModDefModDef[conns[:,4]] .= mod #.= is the broadcasting assignment operator
        while true
            next = next_genes[1]
            if next in genes_remaining
                deleteat!(genes_remaining,findfirst(isequal(next),genes_remaining))
                new_gene = next
                popfirst!(next_genes)
                found = 1
                break
            else if length(next_genes == 1)
                found = 0
                break
            else
                popfirst!(next_genes)
            end
        end
        if found == 0
            break
        end            
    end
end

sOut = [1:length(ModDef) ModDef]

writedlm(outputfile,sOut)
