#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

prefix = args[1]

load("../data/DEtests.RData")
beeSub <- read.csv("../out/substitutions.csv")
beeMK = beeSub[,c(2,6,3,7,4,6,5,7)]
beeMK = cbind(beeMK,rep(22,nrow(beeMK)),rep(1,nrow(beeMK)))

antEvol <- read.csv("../data/MpharAnn.csv")
antEvol = antEvol[!is.na(antEvol$Fixed.Non.Synonymous),]

antSub = antEvol[,c(1,13,17,12,16,15,17,14,16)]

##44 chromosomes sampled (22 workers)
antMK = cbind(antSub[,-c(1)],rep(44,nrow(antSub)),rep(1,nrow(antSub)))

addClass <- function(DEdat,sub,mk,species){
  mk$V10 = rep(1,nrow(mk))
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="queen"]] = 2
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="worker"]] = 3
  write.table(mk,file = paste(prefix,names(DEdat)[2],".",species,".csv",sep = ""),col.names = FALSE,row.names = FALSE,sep=",")
}

sapply(c(2:ncol(beeRes[[2]])),function(x) addClass(beeRes[[2]][,c(1,x)],beeSub,beeMK,"bee"))
sapply(c(2:ncol(antRes[[2]])),function(x) addClass(antRes[[2]][,c(1,x)],antSub,antMK,"ant"))

#Make input from old data 
antMK$class = 1
antMK$class[antEvol$Gaster=="Reproductive"]=2
antMK$class[antEvol$Gaster=="Worker"]=3

write.table(antMK,file = "../MK_alpha_input/old_abd.csv",sep=",",row.names=FALSE,col.names=FALSE)