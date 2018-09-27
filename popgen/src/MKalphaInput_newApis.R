#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

prefix = args[1]
subFile = args[2]
mapFile = args[3]

load("../data/DEtests.RData")
beeSub <- read.csv(subFile)
map = read.table(mapFile,head=F)
colnames(map) = c("Gene","old_gene")
beeSub = merge(beeSub,by="Gene")

mk = beeSub[,c(2,6,3,7,4,6,5,7)]
mk = cbind(mk,rep(496,nrow(mk)),rep(1,nrow(mk)))

addClass <- function(DEdat){
  mk$V10 = rep(1,nrow(mk))
  mk$V10[beeSub$old_gene %in% DEdat$Gene[DEdat[,2]=="queen"]] = 2
  mk$V10[beeSub$old_gene %in% DEdat$Gene[DEdat[,2]=="worker"]] = 3
  write.table(mk,file = paste(prefix,names(DEdat)[2],".bee.csv",sep = ""),col.names = FALSE,row.names = FALSE,sep=",")
}

sapply(c(2:ncol(beeRes[[2]])),function(x) addClass(beeRes[[2]][,c(1,x)],beeSub,beeMK,"bee"))
