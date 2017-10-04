EdgeR <- function(data,design,coef){
  ####Inital steps of standard edgeR analysis
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  ###Calculate FDR
  out <- topTags(diff,n=Inf,adjust.method="BH")$table
  return(out)
}

genFactor <- function(counts){
  factors <- data.frame(sample=colnames(counts))
  factors$stage = 8
  factors$stage[grepl("_E",factors$sample)]=1
  factors$stage[grepl("L1",factors$sample)]=2
  factors$stage[grepl("L2",factors$sample)]=3
  factors$stage[grepl("L3",factors$sample)]=4
  factors$stage[grepl("L4",factors$sample)]=5
  factors$stage[grepl("L5",factors$sample)]=6
  factors$stage[grepl("P",factors$sample)]=7
  factors$tissue="larva"
  factors$tissue[grepl("P",factors$sample)]="pupa"
  factors$tissue[grepl("_E",factors$sample)]="egg"
  factors$tissue[grepl("G\\.",factors$sample)]="gaster"
  factors$tissue[grepl("H\\.",factors$sample)]="head"
  factors$tissue[grepl("M\\.",factors$sample)]="mesosoma"
  factors$NF=NA
  factors$NF[grepl("_N",factors$sample)]="nurse"
  factors$NF[grepl("_F",factors$sample)]="forager"
  factors$caste="worker"
  factors$caste[grepl("_S|_V|_AQ|_G",factors$sample)]="queen"
  factors$caste[grepl("_M",factors$sample)]="male"
  factors$VM=NA
  factors$VM[grepl("_V",factors$sample)]="virgin"
  factors$VM[grepl("_AQ",factors$sample)]="mated"
  factors$colony=1
  factors$colony[grepl(".2",factors$sample)]=2
  factors$colony[grepl(".3",factors$sample)]=3
  for (i in 2:7){
    factors[,i]=as.factor(factors[,i])
  }
  rownames(factors)=factors$sample
  return(factors)
}

tissues <- c("head","gaster","mesosoma")

bee <- read.table("~/GitHub/devnetwork/data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("~/GitHub/devnetwork/data/ants.counts_edit.txt",header=TRUE)
rownames(bee) = bee$gene
rownames(ant) = ant$gene
bee = bee[!grepl("ERCC",bee$gene),-c(1)]
ant = ant[!grepl("ERCC",ant$gene),-c(1)]


ogg2 <- read.csv("~/GitHub/devnetwork/data/HymOGG_hym.csv",sep=" ")
ogg3 <- read.csv("~/GitHub/devnetwork/data/ThreeWayOGGMap.csv",header=TRUE)

factorB <- genFactor(bee)
factorA <- genFactor(ant)

casteDev <- function(caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  design <- model.matrix(~stage+colony,data=droplevels(factors[factors$sample %in% colnames(counts),]))
  out <- EdgeR(counts,design,2:6)
  return(rownames(out)[out$FDR < 0.05])
}

genDevTool <- function(factors,data){
  fQueen <- factors[factors$stage==1|factors$stage==2,]
  fQueen$caste="queen"
  fQueen$sample=paste(fQueen$sample,"_QueenCopy",sep="")
  rownames(fQueen)=fQueen$sample
  factors <- rbind(factors,fQueen) #Add copies of egg and L1 samples for dev toolkit definition
  dataQ <- data[,grepl("L1|_E",colnames(data))]
  colnames(dataQ)=paste(colnames(dataQ),"_QueenCopy",sep="")
  data <- cbind(data,dataQ)#Add copies of egg and L1 samples for dev toolkit definition
  QGenes <- casteDev("queen",data,factors)
  WGenes <- casteDev("worker",data,factors)
  return(QGenes[QGenes %in% WGenes]) #return genes that are differentially expressed across stages in queens and workers
}

BeeDev <- genDevTool(factorB,bee)
AntDev <- genDevTool(factorA,ant)

BeeDevOGG <- ogg2$OGG[ogg2$gene_Amel %in% BeeDev]
AntDevOGG <- ogg2$OGG[ogg2$gene_Mphar %in% AntDev]
TwoSpecDev = BeeDevOGG[BeeDevOGG %in% AntDevOGG]

DevList <- list(BeeDev,AntDev,TwoSpecDev)


##Caste Toolkit
tissueCaste <- function(factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste!="male"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$caste = sample(f$caste,length(f$caste),replace=F)
  }
  design <- model.matrix(~caste+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < 0.05])
}

#############
##Social toolkit
tissueSocial <- function(factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste=="worker"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$NF = sample(f$NF,length(f$NF),replace=F)
  }
  design <- model.matrix(~NF+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < 0.05])
}

save(tissueSocial,tissueCaste,ant,bee,factorA,factorB,DevList,tissues,EdgeR,file="initialvariables.RData")