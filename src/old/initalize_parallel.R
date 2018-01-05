library(edgeR)
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

filterLowly <- function(counts,cut){ #filter out based on counts per million reads
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= ncol(counts)/2
  return(counts[keep,])
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


casteDev <- function(caste,data,factors,fdr){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  design <- model.matrix(~stage+colony,data=droplevels(factors[factors$sample %in% colnames(counts),]))
  out <- EdgeR(counts,design,2:6)
  return(rownames(out)[out$FDR < fdr])
}

genDevTool <- function(factors,data,fdr){
  fQueen <- factors[factors$stage==1|factors$stage==2,]
  fQueen$caste="queen"
  fQueen$sample=paste(fQueen$sample,"_QueenCopy",sep="")
  rownames(fQueen)=fQueen$sample
  factors <- rbind(factors,fQueen) #Add copies of egg and L1 samples for dev toolkit definition
  dataQ <- data[,grepl("L1|_E",colnames(data))]
  colnames(dataQ)=paste(colnames(dataQ),"_QueenCopy",sep="")
  data <- cbind(data,dataQ)#Add copies of egg and L1 samples for dev toolkit definition
  QGenes <- casteDev("queen",data,factors,fdr)
  WGenes <- casteDev("worker",data,factors,fdr)
  return(QGenes[QGenes %in% WGenes]) #return genes that are differentially expressed across stages in queens and workers
}

tissues <- c("head","gaster","mesosoma")

bee <- read.table("~/GitHub/devnetwork/data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("~/GitHub/devnetwork/data/ants.counts_edit.txt",header=TRUE)
rownames(bee) = bee$gene
rownames(ant) = ant$gene
bee = bee[!grepl("ERCC",bee$gene),-c(1)]
ant = ant[!grepl("ERCC",ant$gene),-c(1)]

#Filter out lowly expressed genes
ant <- filterLowly(ant,1)
bee <- filterLowly(bee,1)

ogg2 <- read.csv("~/GitHub/devnetwork/data/HymOGG_hym.csv",sep=" ")
ogg3 <- read.csv("~/GitHub/devnetwork/data/ThreeWayOGGMap.csv",header=TRUE)

factorB <- genFactor(bee)
factorA <- genFactor(ant)

DevWorkflow <- function(fdr){
  
  BeeDev <- genDevTool(factorB,bee,fdr)
  AntDev <- genDevTool(factorA,ant,fdr)
  
  BeeDevOGG <- ogg2$OGG[ogg2$gene_Amel %in% BeeDev]
  AntDevOGG <- ogg2$OGG[ogg2$gene_Mphar %in% AntDev]
  TwoSpecDev = BeeDevOGG[BeeDevOGG %in% AntDevOGG]
  
  DevList <- list(BeeDev,AntDev,TwoSpecDev)
  
  save(tissueSocial,tissueCaste,ant,bee,factorA,factorB,DevList,tissues,EdgeR,ogg2,file=paste("initialvariables",fdr,".RData",sep=""))
  results <- matrix(nrow=1,ncol=9)
  colnames(results) = apply(expand.grid(tissues,c("bee","ant","overlap")), 1, paste, collapse=".")
  write.csv(results,file=paste("casteBoot",fdr,".csv",sep=""))
  write.csv(results,file=paste("socialBoot",fdr,".csv",sep=""))
  return(0)
}

##########
##Caste toolkit
tissueCaste <- function(factors,data,tissue,fdr,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste!="male"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$caste = sample(f$caste,length(f$caste),replace=F)
  }
  design <- model.matrix(~caste+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < fdr])
}

#############
##Social toolkit
tissueSocial <- function(factors,data,tissue,fdr,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste=="worker"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$NF = sample(f$NF,length(f$NF),replace=F)
  }
  design <- model.matrix(~NF+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < fdr])
}

DevWorkflow(0.05)
DevWorkflow(0.1)
DevWorkflow(0.3)
