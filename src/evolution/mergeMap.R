#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

map1 = read.table(args[1],header=F)
map2 = read.table(args[2],header=F)
blast1 = read.table(args[3],header=F)
blast2 = read.table(args[4],header=F)
l1 = read.table(args[5],header=T)
l2 = read.table(args[6],header=T)
t1 = read.table(args[7],header=F)
t2 = read.table(args[8],header=F)

colnames(map1) = colnames(map2) = c("Prot","OGG")
colnames(blast1) = colnames(blast2) = c("transcript","Prot")
names(t1) = names(t2) = c("CDS","transcript")
s1 = merge(map1,blast1,by="Prot")
s2 = merge(map2,blast2,by="Prot")
s1 = merge(s1,t1,by="transcript")
s2 = merge(s2,t2,by="transcript")
s1 = merge(s1,l1,by="CDS")
s2 = merge(s2,l2,by="CDS")

genes = as.list(unique(s1$Gene))
filtered1 <- ldply(lapply(genes,function(x){
  m = s1[s1$Gene==x[1],]
  as.data.frame(m[order(m$Length,decreasing = TRUE),][1,])
}))

genes = as.list(unique(s2$Gene))
filtered2 <- ldply(lapply(genes,function(x){
  m = s2[s2$Gene==x[1],]
  as.data.frame(m[order(m$Length,decreasing = TRUE),][1,])
}))


ogg = merge(filtered1,filtered2,by="OGG")
t = table(ogg$OGG)
keep = names(t)[t==1]
ogg = ogg[ogg$OGG %in% keep,]
t = table(ogg$Gene.x)
keep = names(t)[t==1]
ogg = ogg[ogg$Gene.x %in% keep,]
t = table(ogg$Gene.y)
keep = names(t)[t==1]
ogg = ogg[ogg$Gene.y %in% keep,]



write.table(ogg,file=args[9],row.names=F,col.names=F,quote=F)