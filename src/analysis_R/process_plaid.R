setwd("~/GitHub/devnetwork/")
load("results/PlaidResults.RData")
load("results/collectedPhylo.RData")
library(magrittr)


modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}

beeT <- read.table("data/bees.tpm.txt",header=TRUE)
antT <- read.table("data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

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
  factors$caste = factor(factors$caste,levels = c("queen","male","worker")) #Queen genes will always be down-regulated
  factors$tissue = factor(factors$tissue,levels = c("egg","larva","pupa","head","mesosoma","gaster"))
  factors$NF = factor(factors$NF,levels = c("nurse","forager")) #Make nurse genes down-regulated because nurses should look more like queens (under RGPH)
  return(factors)
}

#Make factors
factorA <- genFactor(ant)
factorB <- genFactor(bee)


freq_set <- function(data,maxMod,factor,tpm,type){
  res <- lapply(data,function(x){
    checkMod = seq(1:nrow(x@NumberxCol))[rowSums(x@NumberxCol) <= maxMod]
    r = list()
    modI = 0
    for (mod in checkMod){
      qp = x@NumberxCol
      samp = factor$sample[x@NumberxCol[mod,]]
      if (sum(grepl(type,samp))==3){
        modI = modI+1
        genes = rownames(tpm)[x@RowxNumber[,mod]]
        r[[modI]]=list(samples=samp,genes=genes)
      }
    }
    if (length(r) > 1){
      rnew = r[[1]]
      rnew$genes = unique(c(r[[2]]$genes,r[[1]]$genes)) #combine genes from queen abdomen module
      r = rnew
    }
    return(r)
  })
  return(res)
}

#Calculate number of significant biclusters for a given set
num_set <- function(set){
  nums = lapply(names(set),function(x){
    length(set[[x]][lapply(set[[x]],length) > 0])
  })
  numDF = data.frame(Type = names(set),Freq=unlist(nums))
}

commonGenes <- function(QGres,tpm){
  df = data.frame(Gene=rownames(tpm),KeepNum=0)
  num = 0
  for (res in QGres){
    if (length(res) == 1){
      num = num+1
      df$KeepNum[df$Gene %in% res[[1]]$genes] = df$KeepNum[df$Gene %in% res[[1]]$genes] + 1
    }
  }
  df$KeepNum = df$KeepNum/num
  return(df)
}

samp_types = sapply(factorA$sample,function(x) gsub(".*_","",x) %>% gsub("\\.[0-9]","",.)) %>%
  unique()

#Output significant biclusters for each sample type
type_bicAmel <- lapply(samp_types,function(x) freq_set(beePl,6,factorB,beeT,x))
type_bicMphar <- lapply(samp_types,function(x) freq_set(antPl,6,factorA,antT,x))
names(type_bicAmel) = names(type_bicMphar) = samp_types
typeNum_Amel = num_set(type_bicAmel)
typeNum_Mphar = num_set(type_bicMphar)


antQG <- type_bicMphar[["AQG"]]
beeQG <- type_bicAmel[["AQG"]]
antQG = antQG[lapply(antQG,length) > 0]
beeQG = beeQG[lapply(beeQG,length) > 0]

aG <- commonGenes(antQG,antT)
bG <- commonGenes(beeQG,beeT)
antG = aG$Gene[aG$KeepNum > 0.9] #Keep genes occurring greater than 60% of the time. Note that there appear to be two peaks in the bee data
beeG = bG$Gene[bG$KeepNum > 0.6]

getConn <- function(tpm,genes,DEtests,DEres){
  tpm = tpm[rownames(tpm) %in% genes,]
  adj1 = cor(t(tpm))^6
  conns = data.frame(Gene=rownames(tpm),conn=rowSums(adj1,na.rm = T))
  bDG =DEres[[1]]
  bDG$DEcat = DEres[[2]]$abdomen
  bFDR = data.frame(Gene=rownames(DEtests[["abdomen"]]),FDR = DEtests[["abdomen"]]$FDR)
  bDG = merge(bDG,bFDR,by="Gene")
  beeGr = bDG[bDG$Gene %in% genes,]
  beeGrC = merge(beeGr,conns,by="Gene")
  beeGrC$DEcat = factor(beeGrC$DEcat,levels = c("queen","worker","nonDE"))
  levels(beeGrC$DEcat) = c("queen","worker","non-DE")
  return(beeGrC)
}

beeGrC = getConn(beeT,beeG,beeTests,beeRes)
antGrC = getConn(antT,antG,antTests,antRes)

write.csv(beeGrC,file = "results/beePlaidGenes.csv")
write.csv(antGrC,file = "results/antPlaidGenes.csv")


