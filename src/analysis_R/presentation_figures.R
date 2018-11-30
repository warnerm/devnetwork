
setwd("~/GitHub/devnetwork/")
load("results/DEtests.RData")
load("results/collectedPhylo.RData")
antConn <- read.csv("results/antConnectivity.csv")
beeConn <- read.csv("results/beeConnectivity.csv")
sexGenes <- read.csv("results/dmel_sexGenes.csv")
library(cowplot)
library(RColorBrewer)
library(VennDiagram)
library(plyr)
library(data.table)
library(zoo)
library(reshape2)
library(EnvStats)
library(grid)
library(magrittr)
library(gridExtra)
library(extrafont)
mypalette2 <- brewer.pal(6,"Blues")
DE_palette = brewer.pal(9,"Blues")
mypalette <- brewer.pal(6,"OrRd")



SexPal = c("firebrick2","slateblue4","gray94")


main_theme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1, color = "black",fill = NA),
        text=element_text(family='Arial'),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 22))

plot2theme <- main_theme + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
                                 axis.title.x = element_text(margin = margin(t=10,b=0,l=0,r=0)),
                                 axis.title.y = element_text(margin = margin(t=0,b=0,l=0,r=10)))



genFactor <- function(counts){
  factors <- data.frame(sample=colnames(counts))
  factors$stage = 8
  factors$stage[grepl("_E",factors$sample)]=1
  factors$stage[grepl("L1",factors$sample)]=2
  factors$stage[grepl("L2",factors$sample)]=3
  factors$stage[grepl("L3",factors$sample)]=4
  factors$stage[grepl("L4",factors$sample)]=5
  factors$stage[grepl("L5",factors$sample)]=6
  factors$stage[grepl("P",factors$sample)]=7
  factors$tissue="larva"
  factors$tissue[grepl("P",factors$sample)]="pupa"
  factors$tissue[grepl("_E",factors$sample)]="egg"
  factors$tissue[grepl("G\\.",factors$sample)]="gaster"
  factors$tissue[grepl("H\\.",factors$sample)]="head"
  factors$tissue[grepl("M\\.",factors$sample)]="mesosoma"
  factors$NF=NA
  factors$NF[grepl("_N",factors$sample)]="nurse"
  factors$NF[grepl("_F",factors$sample)]="forager"
  factors$caste="worker"
  factors$caste[grepl("_S|_V|_AQ|_G",factors$sample)]="queen"
  factors$caste[grepl("_M",factors$sample)]="male"
  factors$VM=NA
  factors$VM[grepl("_V",factors$sample)]="virgin"
  factors$VM[grepl("_AQ",factors$sample)]="mated"
  factors$colony=1
  factors$colony[grepl(".2",factors$sample)]=2
  factors$colony[grepl(".3",factors$sample)]=3
  for (i in 2:7){
    factors[,i]=as.factor(factors[,i])
  }
  rownames(factors)=factors$sample
  factors$caste = factor(factors$caste,levels = c("queen","male","worker")) #Queen genes will always be down-regulated
  factors$tissue = factor(factors$tissue,levels = c("egg","larva","pupa","head","mesosoma","gaster"))
  factors$NF = factor(factors$NF,levels = c("nurse","forager")) #Make nurse genes down-regulated because nurses should look more like queens (under RGPH)
  return(factors)
}


#Take text file and get data frame we want for analysis
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}


#filter out based on counts per million reads
filterLowly <- function(counts,cut,factor){ 
  factor$stc = as.factor(apply(factor[,c(2,3,5)],1,paste,collapse="_"))
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= (ncol(counts)/2)
  keep2 <- lapply(levels(factor$stc),function(x){
    apply(counts[,colnames(counts) %in% factor$sample[factor$stc==x]],1,function(y){
      sum(y > 1)/length(y)
    })
  })
  keepA = do.call(cbind,keep2)
  keepA2 = apply(keepA,1,function(x) sum(x == 1) > 0)
  return(counts[keep|keepA2,])
}

beeT <- read.table("data/bees.tpm.txt",header=TRUE)
antT <- read.table("data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

bee <- read.table("data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("data/ants.counts_edit.txt",header=TRUE)
bee = modifyDF(bee)
ant = modifyDF(ant)

factorA <- genFactor(ant)
factorB <- genFactor(bee)

#############

#Figure 1

#Compare DEG tests across ants and bees
compDef <- function(antR,beeR){
  
  #Define whether or not orthologs exist
  antR$ortholog_found = antR$OGG_found = FALSE
  antR$ortholog_found[antR$Gene %in% AllPS$Gene.x] = TRUE
  antR$OGG_found[antR$Gene %in% ACUogg$gene_Mphar] = TRUE
  aM = melt(antR,id.vars = c("Gene","ortholog_found","OGG_found"))
  aD = ddply(aM,~variable,summarize,
             NDE = sum(value=="nonDE"),
             no_ortholog = sum(value!="nonDE" & !OGG_found),
             OGG = sum(value!="nonDE" & OGG_found))
  
  #Do same thing for apis
  beeR$ortholog_found = beeR$OGG_found = FALSE
  beeR$ortholog_found[beeR$Gene %in% AllPS$Gene.y] = TRUE
  beeR$OGG_found[beeR$Gene %in% ACUogg$gene_Amel] = TRUE
  bM = melt(beeR,id.vars = c("Gene","ortholog_found","OGG_found"))
  bD = ddply(bM,~variable,summarize,
             NDE = sum(value=="nonDE"),
             no_ortholog = sum(value!="nonDE" & !OGG_found),
             OGG = sum(value!="nonDE" & OGG_found))
  colnames(bM)[5] = "value_apis"
  
  #Getting all results together, tabulating
  aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
  bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
  allM = merge(aM,bM,by=c("OGGacu","variable"))
  allD = ddply(allM,~variable,summarize,
               DEboth = sum(value_apis!="nonDE" & value != "nonDE"))
  
  #Calculate number of genes which are DE, have ortholog, and aren't commonly DEG
  aD$DEboth = bD$DEboth = allD$DEboth
  aD$OGG = aD$OGG - aD$DEboth
  bD$OGG = bD$OGG - bD$DEboth
  
  aDM = melt(aD,id.vars = "variable")
  bDM = melt(bD,id.vars = "variable")
  colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")
  aDM$species = "ant"
  bDM$species = "bee"
  
  #Get data back together
  d = rbind(aDM,bDM)
  d$species=as.factor(d$species)
  levels(d$species) = c("M. pharaonis","A. mellifera")
  levels(d$DEtype) = c("NDE"," no ortholog   "," not shared caste-bias   "," shared caste-bias   ")
  return(d)
}

d = compDef(antRes[[2]],beeRes[[2]])

levels(d$stage)[1] = "larva"

levels(d$species) = c("ant","honey bee")
p1m <- ggplot(d[d$DEtype!="NDE",],
              aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",color="black")+
  plot2theme+
  ylim(0,5500)+
  facet_grid(. ~ species)+
  xlab("stage/tissue")+
  scale_fill_manual(values = DE_palette[c(3,7,9)])+
  ylab("number of\ncaste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(panel.spacing = unit(2, "lines"))

#p2 and p3 are plots for the proportions of DGEs
p2 <- ggplot(d[d$DEtype!="NDE" & d$species=="ant",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill",color="black")+
  main_theme+
  scale_fill_manual(values = DE_palette[c(3,7,9)])+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  ylab("proportion")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))


p3 <- ggplot(d[d$DEtype!="NDE" & d$species=="honey bee",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill",color="black")+
  main_theme+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  scale_fill_manual(values = DE_palette[c(3,7,9)])+
  ylab("proportion")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))

#Make proportion plots as insets
pQ <- ggdraw()+
  draw_plot(p1m+
              theme(legend.text = element_text(size=15),
                    legend.key.width = unit(1,"cm"),
                    legend.position = "top"))+
  draw_plot(p2,x=0.18,y=0.53,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.53,height=0.18,width=0.18)

ggsave(pQ,file = "figures/Fig1A.png",height=5,width=9,dpi=300)


levels(d$species) = c("ant","honey bee")
p1m <- ggplot(d[d$DEtype!="NDE",],
              aes(x = stage, y = value,fill=DEtype))+
  geom_bar(stat="identity")+
  plot2theme+
  ylim(0,5500)+
  facet_grid(. ~ species)+
  xlab("stage/tissue")+
  scale_fill_manual(values = rep("black",3))+
  ylab("number of\ncaste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(panel.spacing = unit(2, "lines"))

pQ = p1m+
  theme(legend.text = element_text(size=15),
        legend.key.width = unit(1,"cm"),
        legend.position = "top")
ggsave(pQ,file = "figures/Fig1A_nocolor.png",height=5,width=9,dpi=300)

