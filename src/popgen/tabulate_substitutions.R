#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ann = read.csv(args[1],header=F)
snps = read.csv(args[2],header=F)

library(plyr)

ann = ann[!duplicated(ann),]
snps = snps[!duplicated(snps),]
df <- merge(ann,snps,by="V2")
summed <- ddply(df,~V4,summarize,
                FN = sum(V3.x=="N" & V3.y == "F"),
                FS = sum(V3.x=="S" & V3.y == "F"),
                PN = sum(V3.x=="N" & V3.y == "P"),
                PS = sum(V3.x=="S" & V3.y == "P"))
colnames(summed)[1] = "Gene"
write.csv(summed,file=args[3])