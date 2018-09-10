#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

ann = read.csv(args[1],header=F)
snps = read.csv(args[2],header=F)
silR = read.csv(args[3],header=T)

library(plyr)

ann = ann[!duplicated(ann),]
colnames(ann) = c("CHROM","POS","syn","Gene","transcript")
snps = snps[!duplicated(snps),]
colnames(snps) = c("CHROM","POS","fixed")
df <- merge(ann,snps,by=c("CHROM","POS"))
df = merge(df,silR,by.x="Gene",by.y="gene")
summed <- ddply(df,~Gene,summarize,
                FN = sum(syn=="N" & fixed == "F"),
                FS = sum(syn=="S" & fixed == "F"),
                PN = sum(syn=="N" & fixed == "P"),
                PS = sum(syn=="S" & fixed == "P"),
                Trepl=sum(Trepl),
                Tsil=sum(Tsil))

write.csv(summed,file=args[4],row.names = FALSE)
