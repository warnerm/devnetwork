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

filterLowly <- function(counts,tpm,cut){
  keep = rownames(tpm)[rowSums(tpm > cut) > ncol(tpm)/2]
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

##bootstrapping the whole process
bootOverlap <- function(antF,beeF,antC,beeC,test){
  bee <- list()
  ant <- list()
  TwoSpec <- list()
  
  for (tissue in tissues){
    if (test == "caste"){
      bee[[tissue]]=tissueCaste(beeF,beeC,tissue,scramble=TRUE)
      ant[[tissue]]=tissueCaste(antF,antC,tissue,scramble=TRUE)
    } else {
      bee[[tissue]]=tissueSocial(beeF,beeC,tissue,scramble=TRUE)
      ant[[tissue]]=tissueSocial(antF,antC,tissue,scramble=TRUE)
    }
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% bee[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% ant[[tissue]]]
    TwoSpec[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  AllTest <- list(bee,ant,TwoSpec)
  return(calcOverlap(AllTest))
}

#Calculate overlap of tests
calcOverlap <- function(AllTest){
  overlap = matrix(nrow=1,ncol=9)
  for (i in 1:3){
    for (j in 1:3){
      overlap[1,(i-1)*3+j] = sum(AllTest[[i]][[tissues[j]]] %in% DevList[[i]])
    }
  }
  colnames(overlap) = apply(expand.grid(tissues,c("bee","ant","overlap")), 1, paste, collapse=".")
  return(overlap)
}

#Function to parallelize
runBoots <- function(run){
  while (TRUE){ 
    x <- try(bootOverlap(antF,beeF,antC,beeC,test)) ##Protects from annoying errors that happen with resampling
    if (!inherits(x,"try-error")){
      return(x)
    } 
  }
}

#parallel wrapper
parallelBoot <- function(antF,beeF,antC,beeC,test,nBoots){
  # Initiate cluster; this only works on linux
  cl <- makePSOCKcluster(40,
                         master=system("hostname -i", intern=TRUE))
  
  clusterExport(cl = cl, c(
    "antF","beeF","antC","beeC","test","tissues")) ##Must export these variables for parLapply to see them
  
  # In parallel, go through all permutations
  p <- parLapply(cl,1:nBoots, function(run) {
    runGenie(run)
  })
  stopCluster(cl)
  return(do.call(rbind,p))
}

casteBoot <- parallelBoot(factorA,factorB,ant,bee,"caste",20)
socialBoot <- parallelBoot(factorA,factorB,ant,bee,"social",20)

write.csv(casteBoot,file="casteBoot.csv")
write.csv(socialBoot,file="casteBoot.csv")

