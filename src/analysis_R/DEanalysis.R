######
##Libraries and constants
######

library(edgeR)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(plyr)
library(ppcor)
library(viridis)
library(ggmosaic)
library(lemon)

#Constants
FDR = 0.05
tissues <- c("head","gaster","mesosoma")

######
##Functions
######
#generalized EdgeR workflow, using LR test to identify DEGs
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

#General function to generate fisher test of a continency table
checkOverlap <- function(dat){
  tbl <- rbind(c(dat[1],dat[2]),
               c(dat[3],dat[4]))
  return(fisher.test(tbl))
}

#filter out based on counts per million reads
filterLowly <- function(counts,cut,factor){ 
  factor$stc = as.factor(apply(factor[,c(2,3,5)],1,paste,collapse="_"))
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= (ncol(counts)/2)
  keep2 <- lapply(levels(factor$stc),function(x){
    apply(counts[,colnames(counts) %in% factor$sample[factor$stc==x]],1,function(y){
      sum(y > 1)/length(y)
    })
  })
  keepA = do.call(cbind,keep2)
  keepA2 = apply(keepA,1,function(x) sum(x == 1) > 0)
  return(counts[keep|keepA2,])
}

#generate factors df based on counts fie
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

#Take text file and get data frame we want for analysis
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}

#Add tissue_stage bit
editFactor <- function(f){
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male"&!grepl("_V",f$sample),]) #remove eggs and 1st instars as well as males
  f$tissue_stage = as.factor(do.call(paste,c(f[,c(2,3)],list(sep="_"))))#Want to treat each larval stage separately, as well as each tissue of the adult stage. It's necessary to 
  #do this because the model can't estimate the effect of stage and tissue independently, since there are only adult h/m/g and only larvae of stage 3-6
  relev = levels(f$tissue_stage)
  relev[6:8]=c("8_head","8_mesosoma","8_gaster") #intially alphabetical
  f$tissue_stage = factor(f$tissue_stage,levels=relev)
  return(f)
}

#Treat larva as one stage
editFactor_oneLarv <- function(f){
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male"&!grepl("_V",f$sample),]) #remove eggs and 1st instars as well as males
  f$tissue_stage = as.factor(do.call(paste,c(f[,c(2,3)],list(sep="_"))))#Want to treat each larval stage separately, as well as each tissue of the adult stage. It's necessary to 
  #do this because the model can't estimate the effect of stage and tissue independently, since there are only adult h/m/g and only larvae of stage 3-6
  relev = levels(f$tissue_stage)
  relev[6:8]=c("8_head","8_mesosoma","8_gaster") #intially alphabetical
  f$tissue_stage = factor(f$tissue_stage,levels=relev)
  f$tissue_stage_oneLarv = f$tissue_stage
  f$tissue_stage[f$tissue == "larva"] = "3_larva"
  f$tissue_stage = droplevels(f$tissue_stage)
  return(f)
}


#Return list of nurse/forager DE genes for each tissue
speciesSocial <- function(d,f){
  f = droplevels(f[!is.na(f$NF),])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~NF+colony,data=fs)
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}

#Return list of DE genes at each stage
speciesCaste <- function(d,f,factorFun){
  f = factorFun(f)
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue_stage)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue_stage),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste,data=fs) #bee L4 only came from 1 colony
    results[[lev]] = EdgeR(ds,design,2)
  }
  return(results)
}


#Differential expression based on sex (male vs mated queen)
sexDE <- function(d,f){
  f = droplevels(f[grepl("AQ|_M",f$sample) & f$stage==8,])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste+colony,data=fs) 
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}

#Differential expression based on ovary activation (virgin vs mated queen)
ovaryDE <- function(d,f){
  f = droplevels(f[grepl("AQ|VQ",f$sample),])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~VM+colony,data=fs) 
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}


#Parse social results
parseDE <- function(tests,name1,name2,FDR){
  tests <- lapply(tests,function(x) cbind(x,Gene=rownames(x)))
  all_test <- Reduce(function(x,y) {merge(x,y,by="Gene")},tests)
  all_test2 <- all_test[,c(1,seq(2,(5*length(tests) + 1),by=5))]
  all_test3 <- all_test[,c(1,seq(6,(5*length(tests) + 1),by=5))]
  colnames(all_test2) = c("Gene",names(tests))
  DEtest = all_test2
  DEres <- do.call(cbind,(lapply(1:length(tests),function(j){
    apply(all_test2,1,function(x){
      if (tests[[j]]$FDR[tests[[j]]$Gene==x[1]] < FDR & tests[[j]]$logFC[tests[[j]]$Gene==x[1]] > 0 ) name1
      else if (tests[[j]]$FDR[tests[[j]]$Gene==x[1]] < FDR & tests[[j]]$logFC[tests[[j]]$Gene==x[1]] < 0 ) name2
      else "nonDE"
    })
  })))
  DEtest[,c(2:ncol(all_test2))] = DEres
  return(list(all_test2,DEtest,all_test3))
}

#Identify genes DE across development for the given subset of samples
casteDev <- function(fdr,caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(rownames(out)[out$FDR < fdr])
}


#Generate developmental toolkit: genes DE across develoment in queens and workers
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


############
##Loading in, cleaning data
###########
bee <- read.table("data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("data/ants.counts_edit.txt",header=TRUE)
bee = modifyDF(bee)
ant = modifyDF(ant)

#Make factors
factorA <- genFactor(ant)
factorB <- genFactor(bee)

#Filter out lowly expressed genes
ant <- filterLowly(ant,1,factorA)
bee <- filterLowly(bee,1,factorB)

########
##Identify differentially expressed genes
########
beeTests = speciesCaste(bee,factorB,editFactor)
beeTests_oneLarv = speciesCaste(bee,factorB,editFactor_oneLarv)
antTests = speciesCaste(ant,factorA,editFactor)
antTests_oneLarv = speciesCaste(ant,factorA,editFactor_oneLarv)
ant_sexDE = sexDE(ant,factorA) #Comparing mated queens to males
bee_sexDE = sexDE(bee,factorB)
ant_VM = ovaryDE(ant,factorA) #Comparing mated to virgin queens
bee_VM = ovaryDE(bee,factorB)
beeSocial <- speciesSocial(bee,factorB) #comparing nurses to foragers
antSocial <- speciesSocial(ant,factorA)

#Extract logFC, FDR from DE results
names(antTests_oneLarv) = names(beeTests_oneLarv) = c("larva","pupa","head","thorax","abdomen")
names(beeTests) = names(antTests) = c("L2","L3","L4","L5","pupa","head","thorax","abdomen")

antRes <- parseDE(antTests_oneLarv,"worker","queen",0.1)
beeRes <- parseDE(beeTests_oneLarv,"worker","queen",0.1)
antRes_allstage <- parseDE(antTests,"worker","queen",0.1)
beeRes_allstage <- parseDE(beeTests,"worker","queen",0.1)

#Get social results
names(antSocial) = names(beeSocial) = c("head","thorax","abdomen")
antSocRes <- parseDE(antSocial,"forager","nurse",0.1)
beeSocRes <- parseDE(beeSocial,"forager","nurse",0.1)

#Remove pupae (was including them previously)
DevTool2 <- function(factor,data){
  f = factor[factor$stage != 8 &factor$stage!=7,]
  design <- model.matrix(~stage+colony+caste,data=droplevels(f))
  devGenes <- EdgeR(data[,colnames(data) %in% f$sample],design,2:7)
  return(devGenes)
}

antDevel2 = DevTool2(factorA,ant)
beeDevel2 = DevTool2(factorB,bee)

save(antDevel2,beeDevel2,beeTests,beeTests_oneLarv,antRes,beeRes,antRes_allstage,beeRes_allstage,antSocRes,beeSocRes,antTests,antTests_oneLarv,ant_sexDE,bee_sexDE,ant_VM,bee_VM,beeSocial,antSocial,file = "results/DEtests.RData")
