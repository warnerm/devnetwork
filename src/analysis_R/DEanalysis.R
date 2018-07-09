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

main_theme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1, color = "black",fill = NA),
        text=element_text(family='Arial'),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 22))

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
filterLowly <- function(counts,cut){ 
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= ncol(counts)/2
  return(counts[keep,])
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
parseDE <- function(tests,name1,name2){
  tests <- lapply(tests,function(x) cbind(x,Gene=rownames(x)))
  all_test <- Reduce(function(x,y) {merge(x,y,by="Gene")},tests)
  all_test2 <- all_test[,c(1,seq(2,(5*length(tests) + 1),by=5))]
  all_test3 <- all_test[,c(1,seq(6,(5*length(tests) + 1),by=5))]
  colnames(all_test2) = c("Gene",names(tests))
  DEtest = all_test2
  DEres <- do.call(cbind,(lapply(1:length(tests),function(j){
    apply(all_test2,1,function(x){
      if (tests[[j]]$FDR[tests[[j]]$Gene==x[1]] < 0.1 & tests[[j]]$logFC[tests[[j]]$Gene==x[1]] > 0 ) name1
      else if (tests[[j]]$FDR[tests[[j]]$Gene==x[1]] < 0.1 & tests[[j]]$logFC[tests[[j]]$Gene==x[1]] < 0 ) name2
      else "nonDE"
    })
  })))
  DEtest[,c(2:ncol(all_test2))] = DEres
  return(list(all_test2,DEtest,all_test3))
}


############
##Loading in, cleaning data
###########
bee <- read.table("~/GitHub/devnetwork/data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("~/GitHub/devnetwork/data/ants.counts_edit.txt",header=TRUE)
bee = modifyDF(bee)
ant = modifyDF(ant)

#Filter out lowly expressed genes
ant <- filterLowly(ant,1)
bee <- filterLowly(bee,1)

#Make factors
factorA <- genFactor(ant)
factorB <- genFactor(bee)

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
antRes <- parseDE(antTests_oneLarv,"worker","queen")
beeRes <- parseDE(beeTests_oneLarv,"worker","queen")
antRes_allstage <- parseDE(antTests,"worker","queen")
beeRes_allstage <- parseDE(beeTests,"worker","queen")

#Get social results
antSocRes <- parseDE(antSocial,"forager","nurse")
beeSocRes <- parseDE(beeSocial,"forager","nurse")

save(beeTests,beeTests_oneLarv,antRes,beeRes,antRes_allstage,beeRes_allstage,antSocRes,beeSocRes,antTests,antTests_oneLarv,ant_sexDE,bee_sexDE,ant_VM,bee_VM,beeSocial,antSocial,file = "~/GitHub/devnetwork/data/DEtests.RData")
