library(xgboost)
library(dummies)
setwd("~/GitHub/devnetwork/")
load("results/DEtests.RData")
load("results/collectedPhylo.RData")
tau <- read.csv("results/bee_tau.csv")
sexGenes <- read.csv("results/dmel_sexGenes.csv")
antConn <- read.csv("results/antConnectivity.csv")
beeConn <- read.csv("results/beeConnectivity.csv")
beePi <- read.csv("results/apis.gene.pi.csv")
beeSub <- read.csv("results/substitutions.csv")
antEvol <- read.csv("data/MpharAnn.csv")
antConstraint <- read.table("results/MKtest_constraint_ant.csv")
beeSnipre <- read.csv("results/bayesian_results_apis.csv")
alphaResults <- read.csv("results/collectedAlpha.csv")
colnames(antConstraint) = c("Gene","f")
antSub = antEvol[!is.na(antEvol$Fixed.Non.Synonymous),c(1,13,12,15,14,17,16)]
colnames(antSub) = c("Gene","FN","FS","PN","PS","Trepl","Tsil")
antSub = cbind(antSub,f=as.numeric(as.character((t(antConstraint[2,4:(ncol(antConstraint) - 1)])))))
beeConstraint = read.table("results/MKtest_constraint_bee.csv")
beeSub = cbind(beeSub,f=as.numeric(as.character((t(beeConstraint[2,4:(ncol(beeConstraint) - 1)])))))
antGamma = antEvol[!is.na(antEvol$BSnIPRE.gamma),c(1,19,20)]
beeGamma = beeSnipre[!is.na(beeSnipre$BSnIPRE.gamma),c("gene","BSnIPRE.est","BSnIPRE.gamma")]
colnames(beeGamma)[1] = "Gene"
antSub = merge(antSub,antGamma,by="Gene")
beeAnn = read.csv("results/annotation.csv",header=F) 
beeAnn = beeAnn[!duplicated(beeAnn$V5),]
beeGamma = merge(beeGamma,beeAnn,by.x="Gene",by.y="V5")
beeSub = merge(beeSub,beeGamma[,c(2,3,7)],by.x = "Gene",by.y="V4")

beeT <- read.table("data/bees.tpm.txt",header=TRUE)
antT <- read.table("data/ants.tpm.txt",header=TRUE)
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
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
  factors$tissue[grepl("G\\.",factors$sample)]="abdomen"
  factors$tissue[grepl("H\\.",factors$sample)]="head"
  factors$tissue[grepl("M\\.",factors$sample)]="thorax"
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
  factors$tissue = factor(factors$tissue,levels = c("egg","larva","pupa","head","thorax","abdomen"))
  factors$NF = factor(factors$NF,levels = c("nurse","forager")) #Make nurse genes down-regulated because nurses should look more like queens (under RGPH)
  return(factors)
}

beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]
factorA <- genFactor(antT)
factorB <- genFactor(beeT)

#Calculate euclidean distance
euclDist <- function(res){
  cb = apply(res[,-c(1)],1,function(x) sqrt(sum(x^2))/length(x))
  cb_noAbd = apply(res[,-c(1,ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_noAdult = apply(res[,-c(1,(ncol(res) - 2):ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_larva = apply(res[,-c(1,(ncol(res) - 3):ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_adult = apply(res[,-c(1:(ncol(res) - 3))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_abd = apply(as.data.frame(res[,-c(1:(ncol(res) - 1))]),1,function(x) sqrt(sum(x^2))/length(x))
  results = data.frame(Gene = res$Gene,cb=cb,cb_noAbd=cb_noAbd,cb_noAdult=cb_noAdult,cb_larva=cb_larva,cb_adult = cb_adult,cb_abd=cb_abd)
  return(results)
}

stageExpr <- function(tpm,factor){
  factor$st = as.factor(apply(factor[,c(2,3)],1,paste,collapse="_"))
  factor = droplevels(factor[factor$stage!=2 & factor$stage !=1,])
  expr <- lapply(levels(factor$st),function(x){
    data.frame(rowMeans(tpm[,colnames(tpm) %in% factor$sample[factor$st==x & factor$caste!="male"]]))
  })
  all = do.call(cbind,expr)
  colnames(all) = levels(factor$st)
  all$Gene = rownames(tpm)
  return(all)
}

antExpr <- stageExpr(log(antT^2+1),factorA)
beeExpr <- stageExpr(log(beeT^2+1),factorB)

#For this analysis, we include all larval stages

antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]]) 
antSB = euclDist(antSocRes[[1]])
beeSB = euclDist(beeSocRes[[1]])

aCB = merge(antCB,antSub,by="Gene")
aCB = merge(aCB,Aps,by="Gene")
aCB = merge(aCB,antRes_allstage[[1]],by="Gene")
aCB = merge(aCB,antExpr,by="Gene")
aCB = aCB[aCB$FS > 0 & aCB$PN > 0 & aCB$PS > 0,]
aCB$FN.FS = aCB$FN/(aCB$FS)
aE = data.frame(Gene = rownames(antT),expr=rowMeans(antT))
aCB = merge(aCB,aE)

ant_predictF = aCB[,c(2,4,6,14,20:(ncol(aCB)-1))]
ant_predictFN.FS = aCB[,c(2,4,6,20:ncol(aCB))]
d <- dummy.data.frame(ant_predictFN.FS)
p = get.dummy(d,"psName")
ant_predictF = cbind(ant_predictF[,names(ant_predictF) != "psName"],p)
ggplot(ant_predictF,aes(x=`8_abdomen`,y=FN.FS))+
  geom_hex()+
  geom_smooth()

cor.test(ant_predictF$`8_abdomen`,ant_predictF$FN.FS,method="spearman")
cor.test(aCB$expr,aCB$PS)

dnds = read.table("~/Dropbox/monomorium nurses/data/dnds.txt")
aCB = merge(aCB,dnds[dnds$V2!="STOP",],by.x="Gene",by.y="V1")
aCB$dnds = as.numeric(as.character(aCB$V2))
cor.test(aCB$dnds,aCB$expr,method="spearman")

ggplot(aCB[aCB$dnds <25,],aes(x=`8_abdomen`,y=dnds))+
  geom_hex()+
  geom_smooth()

bee <- read.table("data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("data/ants.counts_edit.txt",header=TRUE)
bee = modifyDF(bee)
ant = modifyDF(ant)
aE = data.frame(Gene = rownames(ant),expr=rowMeans(ant))
aCB = merge(aCB,aE,by="Gene")
cor.test(aCB$FN.FS,aCB$`3_larva`,method="spearman")



#Randomly sample 20% of the data
sample = sample.int(n = nrow(ant_predictF), size = floor(.8*nrow(ant_predictF)), replace = F)
train = ant_predictF[sample,]
test = ant_predictF[-sample,]
train_y = train[,'f']
train_x = train[,names(train) != 'f']
test_y = test[,'f']
test_x = test[,names(test) != 'f']

library(randomForest)
rf_model = randomForest(train_x, y = train_y , ntree = 500, importance = TRUE)
importance_dat = rf_model$importance

sorted_predictors = sort(importance_dat[,1], decreasing=TRUE)



