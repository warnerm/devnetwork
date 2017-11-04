library(edgeR)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
EdgeR <- function(data,design,coef){
  ####Inital steps of standard edgeR analysis
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateDisp(data,design)
  fit <- glmQLFit(dat,design)
  diff <- glmQLFTest(fit, coef=coef) 
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
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
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
    beeCaste[[tissue]]=tissueCaste(factorB,bee,tissue,fdr)
    antCaste[[tissue]]=tissueCaste(factorA,ant,tissue,fdr)
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% beeCaste[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% antCaste[[tissue]]]
    TwoSpecCaste[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  Castelist <- list(beeCaste,antCaste,TwoSpecCaste)
  
  beeNF <- list()
  antNF <- list()
  TwoSpecNF <- list()
  
  for (tissue in tissues){
    beeNF[[tissue]]=tissueSocial(factorB,bee,tissue,fdr)
    antNF[[tissue]]=tissueSocial(factorA,ant,tissue,fdr)
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
  
  png(paste("~/GitHub/devnetwork/results/OverlapTable",fdr,"fdr",".png",sep=""),width=4000,height=1000,res=300)
  grid.table(table)
  dev.off()
  
  write.table(table,file=paste("~/GitHub/devnetwork/results/OverlapTable",fdr,".txt",sep=""),row.names=TRUE)
  return(list(test_list,DevList))
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

############
##Loading bootstrapped results
###########
casteBoot <- read.csv("~/GitHub/devnetwork/data/data.processed/casteBoot.csv")
socialBoot <- read.csv("~/GitHub/devnetwork/data/data.processed/socialBoot.csv")

dfBoot <- function(test,fdr,bootRes,trueRes){
  bootRes <- bootRes[-c(1),-c(1)]
  res <- data.frame(type=rep("bee",9),trueDE=rep(length(trueRes[[1]][["head"]]),9),
                    bootDE=NA,
                    bootC95=NA,
                    bootC99=NA)
  res$type = colnames(bootRes)
  for (i in 1:9){
    res$bootDE[i] = round(mean(bootRes[,i],na.rm=TRUE),3)
    res$bootC95[i] = round(quantile(bootRes[,i],0.95,na.rm=TRUE),3)
    res$bootC99[i] = round(quantile(bootRes[,i],0.99,na.rm=TRUE),3)
  }
  for (i in 1:3){
    for (j in 1:3){
      res$trueDE[(i-1)*3+j] = length(trueRes[[i]][[tissues[[j]]]])
    }
  }
  return(res)
}

caste0.05 <- dfBoot("caste",0.05,casteBoot,fdr_0.05[[1]][[1]])
rownames(caste0.05)=caste0.05$type
plotDf <- melt(caste0.05,id.vars=c("type","bootC95","bootC99"))
plotDf$bootC95[plotDf$variable=="trueDE"]=NA
plotDf$bootC99[plotDf$variable=="trueDE"]=NA
ggplot(plotDf,aes(x=type,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymax=bootC99,ymin=value),position=position_dodge())

caste0.05=caste0.05[,-c(1)]
png("~/GitHub/devnetwork/results/OverlapTableBootCaste.png",width=2000,height=1000,res=300)
grid.table(caste0.05)
dev.off()


caste0.05 <- dfBoot("caste",0.05,socialBoot,fdr_0.05[[1]][[2]])
rownames(caste0.05)=caste0.05$type
caste0.05=caste0.05[,-c(1)]
png("~/GitHub/devnetwork/results/OverlapTableBootSocial.png",width=2000,height=1000,res=300)
grid.table(caste0.05)
dev.off()

###########
##Trying analysis by looking at correlation of F values between Dev toolkit and caste/social toolkit
##########
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


casteDev <- function(caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(out)
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
  worker <- casteDev("worker",data,factors)
  queen <- casteDev("queen",data,factors)
  return(list(worker,queen))
}

getOGG <- function(out,merge.col){
  out$gene = rownames(out)
  out <- merge(out,sharedOGG,by.x="gene",by.y=merge.col)
  return(out)
} 

sharedOGG <- ogg2[ogg2$gene_Mphar %in% rownames(ant) & ogg2$gene_Amel %in% rownames(bee),]
#remove duplicates
sharedOGG$num=1
for (i in 1:nrow(sharedOGG)){
  num = nrow(sharedOGG[sharedOGG$OGG==sharedOGG$OGG[i],])
  sharedOGG$num[i] = num
}
sharedOGG <- sharedOGG[sharedOGG$num==1,]

beeDev <- genDevTool(factorB,bee)
antDev <- genDevTool(factorA,ant)
beeDev_OGG <- lapply(beeDev,getOGG,merge.col="gene_Amel")
antDev_OGG <- lapply(antDev,getOGG,merge.col="gene_Mphar")
Dev_OGG <- c(beeDev_OGG,antDev_OGG)

venn = matrix(nrow=nrow(beeDev_OGG[[1]]),ncol=4)
Fmat = matrix(nrow=nrow(beeDev_OGG[[1]]),ncol=4)
for (i in 1:4){
  for (j in 1:nrow(sharedOGG)){
    row = Dev_OGG[[i]][Dev_OGG[[i]]$OGG %in% sharedOGG$OGG[j],]
    Fmat[j,i] = row$LR
    if (row$FDR < 0.05){
      venn[j,i] = 1
    } else {
      venn[j,i] = 0
    }
  }
}

x <- vennCounts(venn)
png("~/GitHub/devnetwork/results/devVenn.png")
vennDiagram(x,names=c("beeW","beeQ","antW","antQ"),main="Overlap of developmental DE genes")
dev.off()

corMat = matrix(nrow=4,ncol=4)
colnames(corMat) = rownames(corMat) = c("beeW","beeQ","antW","antQ")
for (i in 1:4){
  for (j in 1:4){
    corMat[i,j] = round(cor(Fmat[,i],Fmat[,j]),3)
  }
}
mytheme <- gridExtra::ttheme_default(
  core = list(padding=unit(c(1, 1), "mm"))
)
t <- tableGrob(corMat,theme=mytheme)
png("~/GitHub/devnetwork/results/devCor.png")
grid.arrange(t,top=textGrob("pearson correlation of LR for each gene
            where LR is testing the affect of stage on expression.
            (i.e. whether or not a gene changes expression across pre-pupal development)
            P <<< 0.001 in all cases",vjust=3))
dev.off()

bVenn <- matrix(data=0,nrow=nrow(venn),ncol=2)
bVenn[,1][venn[,1]+venn[,2]==2]=1
bVenn[,2][venn[,3]+venn[,4]==2]=1
x <- vennCounts(bVenn)

chisq.test(rbind(c(x[3,3],x[2,3]),c(x[1,3],x[0,3])))
png("~/GitHub/devnetwork/results/devVenn2.png")
vennDiagram(x,names=c("Amel Dev","Mphar Dev"))
dev.off()

##########
##Caste toolkit
tissueCaste <- function(factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste!="male"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$caste = sample(f$caste,length(f$caste),replace=F)
  }
  design <- model.matrix(~caste+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(out)
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
  return(out)
}

casteList <- list(list())
socList <- list(list())
for (tissue in tissues){
  casteList[["bee"]][[tissue]]=tissueCaste(factorB,bee,tissue)
  socList[["bee"]][[tissue]]=tissueSocial(factorB,bee,tissue)
  casteList[["ant"]][[tissue]]=tissueCaste(factorA,ant,tissue)
  socList[["ant"]][[tissue]]=tissueSocial(factorA,ant,tissue)
}

numDEframe <- function(list){
  beeOGG <- lapply(list[["bee"]],getOGG,merge.col="gene_Amel")
  antOGG <- lapply(list[["ant"]],getOGG,merge.col="gene_Mphar")
  fr = matrix(nrow=3,ncol=3)
  corMat = matrix(nrow=3,ncol=4)
  for (i in 1:3){
    over = merge(beeOGG[[tissues[i]]],antOGG[[tissues[i]]],by="OGG")
    fr[i,1]=sum(beeOGG[[tissues[i]]]$FDR < 0.05)
    fr[i,2]=sum(antOGG[[tissues[i]]]$FDR < 0.05)
    fr[i,3]=sum(over$FDR.x < 0.05 & over$FDR.y < 0.05)
    mat = cor.test(over$logFC.x,over$logFC.y)
    corMat[i,1] = signif(mat$estimate,4)
    corMat[i,2] = signif(mat$p.value,4)
    mat = cor.test(abs(over$logFC.x),abs(over$logFC.y))
    corMat[i,3] = signif(mat$estimate,4)
    corMat[i,4] = signif(mat$p.value,4)
  }
  rownames(fr)= rownames(corMat) = tissues
  colnames(fr) = c("bee","ant","overlap")
  colnames(corMat) = c("R","pVal","R_abs","pVal_abs")
  return(list(fr,corMat))
}

cDE <- numDEframe(casteList)
sDE <- numDEframe(socList)

res = c(cDE,sDE)
l = lapply(res,tableGrob)
png("~/GitHub/devnetwork/results/casteTables.png",height=1000,width=3000,res=300)
grid.arrange(l[[1]],l[[3]],l[[2]],l[[4]])
dev.off()

devGenes <- function(dev){
  dev2 <- merge(dev[[1]],dev[[2]],by="OGG")
  return(dev2$OGG[dev2$FDR.x < 0.05 & dev2$FDR.y < 0.05])
}

devLR <- function(dev){
  dev2 <- merge(dev[[1]],dev[[2]],by="OGG")
  dev2$meanLR = (dev2$LR.x+dev2$LR.y)/2
  return(dev2)
}
dev = list(beeDev_OGG,antDev_OGG)
devGene <- lapply(dev,devGenes)
devLRlist <- lapply(dev,devLR)

numDEoverlap <- function(list,test){
  beeOGG <- lapply(list[["bee"]],getOGG,merge.col="gene_Amel")
  antOGG <- lapply(list[["ant"]],getOGG,merge.col="gene_Mphar")
  nDE <- matrix(nrow=3,ncol=9)
  corMat <- matrix(nrow=3,ncol=4)
  for (i in 1:3){
    beeG = beeOGG[[tissues[i]]][beeOGG[[tissues[i]]]$FDR < 0.05,]
    antG = antOGG[[tissues[i]]][antOGG[[tissues[i]]]$FDR < 0.05,]
    overG = beeG$OGG[beeG$OGG %in% antG$OGG]
    nDE[i,1] = length(devGene[[1]])
    nDE[i,2] = nrow(beeG)
    nDE[i,3] = sum(beeG$OGG %in% devGene[[1]])
    nDE[i,4] = length(devGene[[2]])
    nDE[i,5] = nrow(antG)
    nDE[i,6] = sum(antG$OGG %in% devGene[[2]])
    nDE[i,7] = sum(devGene[[1]] %in% devGene[[2]])
    nDE[i,8] = length(overG)
    nDE[i,9] = sum(overG %in% devGene[[1]][devGene[[1]] %in% devGene[[2]]])
    mat = cor.test(devLRlist[[1]]$meanLR,beeOGG[[tissues[i]]]$LR)
    corMat[i,1] = signif(mat$estimate,4)
    corMat[i,2] = signif(mat$p.value,4)
    mat = cor.test(devLRlist[[2]]$meanLR,antOGG[[tissues[i]]]$LR)
    corMat[i,3] = signif(mat$estimate,4)
    corMat[i,4] = signif(mat$p.value,4)
  }
  rownames(nDE) = rownames(corMat) = tissues
  colnames(nDE) =  apply(expand.grid(c("devel",test,paste("devel_",test,sep="")),c("bee","ant","overlap")), 1, paste, collapse=".")
  colnames(corMat) = c("bee_R","bee_pval","ant_R","ant_pval")
  return(list(nDE,corMat))
}

cOver <- numDEoverlap(casteList,"caste")
sOver <- numDEoverlap(socList,"social")

tbl = list(cOver[[1]],sOver[[1]])
tbls <- lapply(tbl,tableGrob)

png("~/GitHub/devnetwork/results/overTables.png",height=1500,width=4000,res=300)
grid.arrange(tbls[[1]],tbls[[2]])
dev.off()

##############
##Comparing to genes defined as developmental toolkit in D. melanogaster
#############
dmel <- read.table("~/GitHub/devnetwork/data/DmelDevelOGG.txt")
oggDev = droplevels(ogg3[ogg3$OGG %in% dmel$V1,])
t = table(oggDev$OGG)
t = t[t==1]
oggDevF = oggDev[oggDev$OGG %in% names(t),] #Filter for 1-1-1 orthologs (developmental genes)

t = table(ogg3$OGG)
t = t[t==1]
ogg111 = ogg3[ogg3$OGG %in% names(t),] #Filter for 1-1-1 orthologs (all genes)
ogg111 = ogg111[ogg111$gene_Mphar %in% rownames(antDev[[1]]) &
                  ogg111$gene_Amel %in% rownames(beeDev[[1]]),] #Filter for oggs in the analysis


beeDW = rownames(beeDev[[1]])[beeDev[[1]]$FDR < 0.05]
antDW = rownames(antDev[[1]])[antDev[[1]]$FDR < 0.05]
beeDQ = rownames(beeDev[[2]])[beeDev[[2]]$FDR < 0.05]
antDQ = rownames(antDev[[2]])[antDev[[2]]$FDR < 0.05]

beeD = beeDW[beeDW %in% beeDQ]
antD = antDW[antDW %in% antDQ]

beeD_ogg = ogg111$OGG[ogg111$gene_Amel %in% beeD]
antD_ogg = ogg111$OGG[ogg111$gene_Mphar %in% antD]
DevOGG = beeD_ogg[beeD_ogg %in% antD_ogg]

MasterDev = DevOGG[DevOGG %in% oggDevF$OGG] #OGGs identified as dev toolkit by transcriptomic and Dmel annotation

a = rbind(c(length(MasterDev),nrow(oggDevF) - length(MasterDev)),
          c(length(DevOGG) - length(MasterDev),
            nrow(ogg111) - nrow(oggDevF) - length(DevOGG) + length(MasterDev)))

fisher.test(a,alternative="greater") #They is less association than expected
oggDevF[oggDevF$OGG %in% MasterDev,]

##############
##Are social and caste toolkits enriched for developmental genes, as defined in Dmel?
##############
getOGG <- function(list, column){
  gene = rownames(list)[list$FDR < 0.05]
  oggC = ogg111$OGG[ogg111[,column] %in% gene]
  return(oggC)
}

checkOver <- function(oggC,numC){
  MasterOGG = oggDevF$OGG[oggDevF$OGG %in% oggC]
  a = rbind(c(length(MasterOGG),nrow(oggDevF) - length(MasterOGG)),
            c(length(oggC) - length(MasterOGG),
              nrow(ogg111) - nrow(oggDevF) - length(oggC) + length(MasterOGG)))
  r = fisher.test(a,alternative="greater")
  return(c(length(oggC),length(MasterOGG),r$p.value))
}

res = list(list(list()))
tests = c("caste","social")
bigL = list(casteList,socList)
names(bigL) = tests
cols = c(2,3)
names(cols) = c("ant","bee")

for (test in tests){
  for (spec in c("bee","ant")){
    for (tissue in tissues){
      oggC = getOGG(bigL[[test]][[spec]][[tissue]],cols[spec])
      res[[test]][[spec]][[tissue]]=checkOver(oggC)
    }  
  }
}









