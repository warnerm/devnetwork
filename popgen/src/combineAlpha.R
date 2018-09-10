#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

inDir = args[1]
outputFile = args[2]

collect <- function(stage,species){
  df <- read.csv(paste(inDir,"alpha_locusF_",stage,".",species,".csv",sep=""),sep="\t")
  df = df[-(nrow(df)),]
  df = df[,c("LnL","theta","tau","alpha.Cat1","alpha.Cat2","alpha.Cat3")]
  df = df[-c(1,2),] 
  #alpha = apply(df[,c(4:6)],2,function(x) mean(as.numeric(as.character(x)))) #First line is the true values (no permutation test)
  #c1 = apply(df[,c(4:6)],2,function(x) quantile(x,0.025,na.rm=T)) #Note that every result is printed twice, but we don't really care because we are just looking for quantiles anyway
  #c2 = apply(df[,c(4:6)],2,function(x) quantile(x,0.975,na.rm=T))
  #d = as.data.frame(do.call(rbind,list(alpha,c1,c2)))
  d = df
  d$stage=stage
  d$species=species
  return(d)
}

res = ldply(lapply(c("ant","bee"), function(x){
  ldply(lapply(c("larva","pupa","head","thorax","abdomen"), function(y){
    collect(y,x)
  }))
}))

write.csv(res,file=outputFile,row.names=FALSE)

