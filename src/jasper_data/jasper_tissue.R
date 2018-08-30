setwd("~/GitHub/devnetwork/src/")

counts <- read.csv("../data/Jasper_counts.csv")
fpkm <- read.csv("../data/Jasper_fpkm.csv")
samples <- read.table("../data/samples_jasper.txt",head=T)
rownames(counts) = counts$gene_id
counts = counts[,-c(1)]
rownames(fpkm) = fpkm$gene_id
fpkm = fpkm[,-c(1)]
fpkm = fpkm[rowSums(fpkm) > 0,]

calculateTau <- function(factor,expr){
  meanExpr <- lapply(levels(factor$tissue),function(x) 
    rowSums(as.data.frame(expr[,colnames(expr) %in% factor$SRA[factor$tissue==x]]))/sum(factor$tissue==x))
  meanExpr <- as.data.frame(do.call(cbind,meanExpr))
  colnames(meanExpr) = levels(factor$sample)
  geneDeviation = apply(meanExpr,1,function(x){
    sum(ldply(lapply(c(1:12),function(i) 1-x[i]/max(x))))
  })
  
  tau = geneDeviation/(length(levels(factor$tissue))-1)
  
  return(list(tau,meanExpr))
}

beeTau = calculateTau(samples,fpkm)

tau <- data.frame(Gene = names(beeTau[[1]]),tau=tau)

write.csv(tau,file = "../results/bee_tau.csv")
# 
# TGmap <- read.table("../phylostratigraphy/out/TGmap_Amel.txt")
# TNmap <- as.data.frame(fread("../data/AmelTranName.txt",sep="~",header=FALSE))
# 
# AmelName <- merge(TGmap,TNmap,by.x = "V2",by.y = "V1")[,c(2,3)]
# colnames(AmelName) = c("Gene","GeneName")
# AmelName$GeneName = gsub(" isoform X[0-9]","",AmelName$GeneName)
# aName = AmelName[!duplicated(AmelName$Gene),]
# aName = merge(aName,DmelSC,by.x="Gene",by.y="gene_Amel",all.x=T)
# aName = aName[,c(1,2,3)]
# colnames(aName) = c("gene_Amel","geneName_Amel","gene_Dmel")
# aName = aName[!duplicated(aName$gene_Amel),]
# write.table(aName,file="~/GitHub/workshop-iussi2018/Apis_mellifera_geneName.txt")
