library(lme4)
library(R2jags)
library(arm)
library(R2WinBUGS)

print("Loading snipre source")
source("snipre/B_SnIPRE_source.R")
source("snipre/my.jags2.R")

print("Loading data")
#load data

annotation <- read.csv('../data/popgen/var/annotation.csv',header=F) # read data on mutation effect
colnames(annotation) <- c("chrom", "pos", "effect", "loc", "rna")
annotation$pos <- paste(annotation$chrom,annotation$pos)
annotation <- annotation[,-1]
snps <- read.csv('../data/popgen/var/snps.csv',header=F) # read fixed vs polymorphic
colnames(snps) <- c("chrom", "pos", "state")
snps$pos <- paste(snps$chrom,snps$pos)
snps <- snps[,-1]
silentReplacement <- read.csv('../data/popgen/var/silentReplacement.csv', header=T) # read expected silent replacement

print("Data loaded")
byPos <- merge(snps,annotation, by.x="pos", by.y="pos")

rownames(silentReplacement) <- silentReplacement$isoform
silentReplacement <- silentReplacement[,-1]

mk <- data.frame(FS = as.data.frame(table(subset(byPos,state=="F" & effect=="S")$rna))$Freq,
FR = as.data.frame(table(subset(byPos,state=="F" & effect=="N")$rna))$Freq,
PS = as.data.frame(table(subset(byPos,state=="P" & effect=="S")$rna))$Freq,
PR = as.data.frame(table(subset(byPos,state=="P" & effect=="N")$rna))$Freq)
rownames(mk) <- as.data.frame(table(subset(byPos,state=="F" & effect=="S")$rna))$Var1

dnds <- merge(mk,silentReplacement,by=0)
colnames(dnds)[1] <- "gene"
dnds$nout <- 2
dnds$npop <- 20*2
#adjust column order
dnds <- dnds[,c("gene", "FS", "FR", "PS", "PR", "Tsil", "Trepl", "nout", "npop")]

print("starting snipre computation")

bugs.directory='~/.wine/drive_c/Program Files (x86)/WinBUGS14'
BSnIPRE.run(dnds, burnin = 10000, thin = 4, iter = 15000)
load("samples")
res <- BSnIPRE(samples, dnds)
bres <- res$new.dataset
write.csv(bres, file = "../out/bayesian_results.csv", row.names = FALSE)


