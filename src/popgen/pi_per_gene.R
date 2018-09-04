#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ann = read.csv(args[1],header=F)
pi_site = read.csv(args[2],header=T)
pi_site = pi_site[!duplicated(pi_site),]
pi_ann = merge(pi_site,ann,by.x=c("CHROM","POS"),by.y=c("V1","V2"))
pi_mean = ddply(pi_site,~ V4,summarize,
                pi = mean(PI))
colnames(pi_mean)[1] = "Gene"

write.csv(pi_mean,file=args[2],row.names = FALSE)