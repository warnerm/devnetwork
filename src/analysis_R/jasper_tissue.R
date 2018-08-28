setwd("~/GitHub/devnetwork/src/")

counts <- read.csv("../data/Jasper_counts.csv")
fpkm <- read.csv("../data/Jasper_fpkm.csv")
samples <- read.table("../data/samples_jasper.txt",head=T)
rownames(counts) = counts$gene_id
counts = counts[,-c(1)]
rownames(fpkm) = fpkm$gene_id
fpkm = fpkm[,-c(1)]
fpkm = fpkm[rowSums(fpkm) > 0,]

calculatetauSoc <- function(factor,expr){
  factor$stc = as.factor(apply(factor[,c(2,3,5)],1,paste,collapse="_"))
  meanExpr <- lapply(levels(factor$stc),function(x) 
    rowSums(as.data.frame(expr[,colnames(expr) %in% factor$sample[factor$stc==x]]))/sum(factor$stc==x))
  meanExpr <- as.data.frame(do.call(cbind,meanExpr))
  colnames(meanExpr) = levels(factor$stc)
  geneDeviation = apply(meanExpr,1,function(x){
    sum(ldply(lapply(c(1:12),function(i) 1-x[i]/max(x))))
  })
  
  tau = geneDeviation/(length(levels(factor$stc))-1)
  
  return(list(tau,meanExpr))
}

antTau = calculatetauSoc(factorA,antT[rowSums(antT) > 0,])
beeTau = calculatetauSoc(factorB,beeT[rowSums(beeT) > 0,])

save(beeTau,antTau,file = "../out/tau_results.RData")

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
