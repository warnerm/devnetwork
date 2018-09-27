#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

prefix = args[1]
outputFile = args[2]

collect <- function(stage,species){
  df <- read.csv(paste(prefix,stage,".",species,".csv",sep=""),sep="\t",head=T)
  #df = df[-c(1:3),]
  #df = df[-seq(1,nrow(df) - 1, by=2),c("alpha.Cat1","alpha.Cat2","alpha.Cat3")]
  df$stage=stage
  df$species=species
  return(as.data.frame(df[1,]))
}

res = ldply(lapply(c("bee"), function(x){
  ldply(lapply(c("larva","pupa","head","thorax","abdomen"), function(y){
    collect(y,x)
  }))
}))

write.csv(res,file=outputFile,row.names=FALSE)

