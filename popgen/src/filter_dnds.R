#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

dnds = read.table(args[1],header=T)
lengths = read.table(args[2],header=T)
map = read.table(args[3],header=F,sep=" ")

d = cbind(dnds,map)
d = merge(d,lengths,by.x="V1",by.y="GeneID")

d = d[d$dN_dS!="STOP" & d$dN_dS!= "SHORT" & d$dN_dS!=99.000,]
d$dN_dS = as.numeric(as.character(d$dN_dS))
d = d[d$dN_dS!=99.000,]

genes = as.list(unique(d$Gene))
filtered <- ldply(lapply(genes,function(x){
  m = d[d$Gene==x[1],]
  as.data.frame(m[order(m$Length,decreasing = TRUE),][1,])
}))

write.table(filtered,file=args[4])