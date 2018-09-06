#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

inDir = args[1]
outputFile = args[2]

collect <- function(stage,species){
  df <- read.csv(paste(inDir,"alpha_globalF_",stage,".",species,".csv",sep=""),sep="\t")
  alpha = apply(df[,c(5:7)],2,mean)
  c1 = apply(df[,c(5:7)],2,function(x) quantile(x,0.025))
  c2 = apply(df[,c(5:7)],2,function(x) quantile(x,0.975))
  d = as.data.frame(do.call(cbind,list(alpha,c1,c2)))
  d$caste = c("nonDE","worker","queen")
  d$stage=stage
  d$species=species
  colnames(d)[1:3] = c("alpha","c1","c2")
  return(d)
}

res = ldply(lapply(c("ant","bee"), function(x){
  ldply(lapply(c("L2","L3","L4","L5","pupa","head","thorax","abdomen"), function(y){
    collect(y,x)
  }))
}))

write.csv(res,file=outputFile,row.names=FALSE)

