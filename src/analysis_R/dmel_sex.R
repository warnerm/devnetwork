#Load in endopterygota definitions
load("../phylostratigraphy/out/collectedPhylo.RData")

#Identify sex-associated genes in Drosophila based on RNA-seq data
counts <- read.csv("../data/counts_drosophila.csv")
fpkm <- read.csv("../data/fpkm_drosophila.csv")

factors <- read.csv("../data/drosophila_sra.csv")
rownames(counts) = rownames(fpkm) = counts$gene_id
counts = counts[,-c(1)]
fpkm = fpkm[,-c(1)]

keep = apply(fpkm,1,function(x) sum(x > 1) >= length(x)/2)
counts = counts[keep,]

factors$time_factor = as.factor(factors$time)

#Sex-bias
design <- model.matrix(~sex+time_factor,data=droplevels(f_sex[f_sex$time>25,]))
d = counts[,colnames(counts) %in% f_sex$SRA[f_sex$time>25]]
sexGenes <- EdgeR(d,design,2)
sexGenes$Gene = rownames(sexGenes)
sexGenes = merge(sexGenes,ENDogg,by.x="Gene",by.y="gene_Dmel")

write.csv(sexGenes,file="../out/dmel_sexGenes.csv")