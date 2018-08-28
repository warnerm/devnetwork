setwd("~/GitHub/devnetwork/src/")
library(edgeR)

EdgeR <- function(data,design,coef){
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  ##Calculates FDR
  out <- topTags(diff,n=Inf,adjust.method="BH")$table
  return(out)
}

#Load in endopterygota definitions
load("../results/collectedPhylo.RData")

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
design <- model.matrix(~sex+time_factor,data=droplevels(factors[factors$time>25,]))
d = counts[,colnames(counts) %in% factors$SRA[factors$time>25]]
sexGenes <- EdgeR(d,design,2)
sexGenes$Gene = rownames(sexGenes)
sexGenes = merge(sexGenes,ENDogg,by.x="Gene",by.y="gene_Dmel")

write.csv(sexGenes,file="../results/dmel_sexGenes.csv")