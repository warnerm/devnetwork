#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(plyr)

ann = read.csv(args[1],header=F)
pi_site = read.table(args[2],header=T)
pi_site = pi_site[!duplicated(pi_site),]
ann = ann[!duplicated(ann),]
colnames(ann) = c("CHROM","POS","TYPE","Gene","Transcript")
pi_ann = merge(pi_site,ann,by=c("CHROM","POS"))
pi_mean = ddply(pi_ann,~Gene,summarize,
                pi = mean(PI))

write.csv(pi_mean,file=args[3],row.names = FALSE)