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

write.csv(factorB,"~/GitHub/devnetwork/data/factors_bees.csv")
write.csv(factorA,"~/GitHub/devnetwork/data/factors_ants.csv")

#########Functions for the analysis
#############
###Devel toolkit
casteDev <- function(fdr,caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  design <- model.matrix(~stage+colony,data=droplevels(factors[factors$sample %in% colnames(counts),]))
  out <- EdgeR(counts,design,2:6)
  return(rownames(out)[out$FDR < fdr])
}

genDevTool <- function(fdr,factors,data){
  fQueen <- factors[factors$stage==1|factors$stage==2,]
  fQueen$caste="queen"
  fQueen$sample=paste(fQueen$sample,"_QueenCopy",sep="")
  rownames(fQueen)=fQueen$sample
  factors <- rbind(factors,fQueen) #Add copies of egg and L1 samples for dev toolkit definition
  dataQ <- data[,grepl("L1|_E",colnames(data))]
  colnames(dataQ)=paste(colnames(dataQ),"_QueenCopy",sep="")
  data <- cbind(data,dataQ)#Add copies of egg and L1 samples for dev toolkit definition
  QGenes <- casteDev(fdr,"queen",data,factors)
  WGenes <- casteDev(fdr,"worker",data,factors)
  return(QGenes[QGenes %in% WGenes]) #return genes that are differentially expressed across stages in queens and workers
}

##########
##Caste toolkit
tissueCaste <- function(fdr,factors,data,tissue,scramble=FALSE){
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
tissueSocial <- function(fdr,factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste=="worker"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$NF = sample(f$NF,length(f$NF),replace=F)
  }
  design <- model.matrix(~NF+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < fdr])
}

##Entire workflow with different fdr values
workflow <- function(fdr){
  BeeDev <- genDevTool(fdr,factorB,bee)
  AntDev <- genDevTool(fdr,factorA,ant)
  
  BeeDevOGG <- ogg2$OGG[ogg2$gene_Amel %in% BeeDev]
  AntDevOGG <- ogg2$OGG[ogg2$gene_Mphar %in% AntDev]
  TwoSpecDev = BeeDevOGG[BeeDevOGG %in% AntDevOGG]
  
  DevList <- list(BeeDev,AntDev,TwoSpecDev)
  nGene <- list(nrow(bee),nrow(ant),length(commonOGG))
  
  beeCaste <- list()
  antCaste <- list()
  TwoSpecCaste <- list()
  
  for (tissue in tissues){
    beeCaste[[tissue]]=tissueCaste(fdr,factorB,bee_filt,tissue)
    antCaste[[tissue]]=tissueCaste(fdr,factorA,ant_filt,tissue)
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% beeCaste[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% antCaste[[tissue]]]
    TwoSpecCaste[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  Castelist <- list(beeCaste,antCaste,TwoSpecCaste)
  
  beeNF <- list()
  antNF <- list()
  TwoSpecNF <- list()
  
  for (tissue in tissues){
    beeNF[[tissue]]=tissueSocial(fdr,factorB,bee,tissue)
    antNF[[tissue]]=tissueSocial(fdr,factorA,ant,tissue)
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% beeNF[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% antNF[[tissue]]]
    TwoSpecNF[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  NFlist <- list(beeNF,antNF,TwoSpecNF)
  
  #######
  ##Hypothesis tests
  ######Make table with # of genes in each of the lists above (not developmental), plus number of genes in analysis (in parentheses)
  ######Then add number of genes in each overlap comparison
  ######Rows are head, mesosoma, and gaster for both caste and social
  ######Can do a quick chi square, but also do bootstrapping for confidence intervals based on scrambled data
  test_list <- list(Castelist,NFlist)
  
  table <- matrix(nrow = 8, ncol = 9)
  for (i in 1:2){
    for (j in 1:3){
      for (k in 1:3){
        table[(i-1)*3+j,k] = length(test_list[[i]][[k]][[tissues[j]]])
      }
    }
  }
  
  for (i in 1:3){
    table[7,i] = length(DevList[[i]])
    table[7,i+3] = length(DevList[[i]])
    table[7,i+6] = length(DevList[[i]])
    table[8,i] = nGene[[i]]
    table[8,i+3] = nGene[[i]]
    table[8,i+6] = nGene[[i]]
  }
  
  
  
  for (i in 1:3){
    for (j in 1:2){
      for (k in 1:3){
        table[(j-1)*3+k,i+3] = sum(DevList[[i]] %in% test_list[[j]][[i]][[tissues[k]]])
        chi = rbind(c(table[(j-1)*3+k,i+3],
                      table[(j-1)*3+k,i]-table[(j-1)*3+k,i+3]),
                    c(table[7,i] - table[(j-1)*3+k,i+3],
                      nGene[[i]] - table[7,i] - table[(j-1)*3+k,i] + table[(j-1)*3+k,i+3]))
        table[(j-1)*3+k,i+6] = signif(fisher.test(chi,alternative="greater")$p.value,digits=3)
      }
    }
  }
  
  rownames(table) <- c(apply(expand.grid(tissues,c("Caste","NF")), 1, paste, collapse="."),"DevelToolkit","nGene_Analysis")
  colnames(table) <- c("Amel","Mpharaonis","TwoSpecies","Amel_Devel",
                       "Mpharaonis_Devel","TwoSpecies_Devel",
                       "Amel_fisherP","Mpharaonis_fisherP","TwoSpecies_fisherP")
  
  png(paste("~/GitHub/devnetwork/results/OverlapTable","fdr",".png",sep=""),width=4000,height=1000,res=300)
  grid.table(table)
  dev.off()
  
  write.table(table,file=paste("~/GitHub/devnetwork/results/OverlapTable",fdr,".txt",sep=""),row.names=TRUE)
  return(list(test_list,DevList))
}

#constants
tissues <- c("head","gaster","mesosoma")

BeeOGG <- ogg2$OGG[ogg2$gene_Amel %in% rownames(bee)]
AntOGG <- ogg2$OGG[ogg2$gene_Mphar %in% rownames(ant)]
commonOGG <- BeeOGG[BeeOGG %in% AntOGG]

#############
##Running program
#############
fdr_0.05 <- workflow(0.05)
fdr_0.1 <- workflow(0.1)
fdr_0.3 <- workflow(0.3)


