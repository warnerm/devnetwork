
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
PS_palette = brewer.pal(9,"Greens")

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

levels(d$stage)[1] = "larva*"

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
  draw_plot(p2,x=0.16,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.57,y=0.58,height=0.18,width=0.18)

#Figure 1b- nurse vs forager
d = compDef(antSocRes[[2]],beeSocRes[[2]])
levels(d$DEtype) = c("NDE"," no ortholog   "," not shared task-bias   "," shared task-bias   ")

levels(d$species) = c("ant","honey bee")
p1m <- ggplot(d[d$DEtype!="NDE",],
              aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",color="black")+
  plot2theme+
  ylim(0,5500)+
  facet_grid(. ~ species)+
  xlab("tissue")+
  scale_fill_manual(values = DE_palette[c(3,7,9)])+
  ylab("number of\ntask-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "none",
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(panel.spacing = unit(2, "lines"))

#p2 and p3 are proportion plots
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

pSoc <- ggdraw()+
  draw_plot(p1m+
              theme(legend.text = element_text(size=15),
                    legend.key.width = unit(1,"cm"),
                    legend.position = "top"))+
  draw_plot(p2,x=0.16,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.57,y=0.58,height=0.18,width=0.18)

#Calculate proportion of each age class in DEG categories, focusing on abdomen
psO = merge(AllPS_sum,ACUogg,by="OGGacu")
psO$a2 = "not caste-biased"
psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen!="non-DE"] | 
       psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen!="non-DE"]] = "inconsistent"

psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="worker"] & 
         psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen=="worker"]] = "shared worker"
psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="queen"] & 
         psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen=="queen"]] = "shared queen"

levels(psO$psName) = c("ancient","insect","hymenopteran","aculeate")

a = rbind(c(sum(grepl("shared",psO$a2) & psO$psName=="ancient"),
            sum(grepl("shared",psO$a2) & psO$psName!="ancient")),
          c(sum(grepl("not caste-biased",psO$a2) & psO$psName=="ancient"),
            sum(grepl("not caste-biased",psO$a2) & psO$psName!="ancient")))
a = rbind(c(sum(grepl("shared queen",psO$a2) & psO$psName=="ancient"),
            sum(grepl("shared queen",psO$a2) & psO$psName!="ancient")),
          c(sum(grepl("shared worker",psO$a2) & psO$psName=="ancient"),
            sum(grepl("shared worker",psO$a2) & psO$psName!="ancient")))

psO$a2 = factor(psO$a2,levels = rev(c("shared queen","shared worker","inconsistent","not caste-biased")))

levels(psO$psName) = paste(" ", levels(psO$psName),"   ",sep="")


p <- ggplot(psO,aes(x = a2,fill = forcats::fct_rev(psName)))+
  ylab("proportion")+
  geom_bar(stat = "count",position = "fill",color="black",alpha=0.8)+
  scale_fill_manual(values = PS_palette[c(2,4,6,8)],name = "phylostrata")+
  main_theme+
  coord_flip()+
  xlab("abdominal DEG type")+
  scale_x_discrete(expand=c(0,0))+
  guides(fill = guide_legend(reverse=T))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle=-25,hjust=0),
        legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5,2,0.5,2),"cm"))

ggsave(arrangeGrob(pQ,pSoc,p,heights = c(0.37,0.37,0.26)),file="figures/Fig1.png",height=14,width=10,dpi=1200)

###################

#Figure 2
antPlaid <- read.csv("results/antPlaidGenes.csv")
beePlaid <- read.csv("results/beePlaidGenes.csv")
antGrC = merge(antPlaid,Aps,by="Gene")
beeGrC = merge(beePlaid,Bps,by="Gene")
antGrC$species= "ant"
beeGrC$species="honey bee"

#Scale connectivity to the top value
antGrC$conn = antGrC$conn/max(antGrC$conn)
beeGrC$conn = beeGrC$conn/max(beeGrC$conn)

#Identify genes that are present in the bicluster of both species
antGrC$conBic=beeGrC$conBic="non-conserved"
antGrC$conBic[antGrC$Gene %in% oB$gene_Mphar]="conserved"
beeGrC$conBic[beeGrC$Gene %in% oA$gene_Amel]="conserved"

#Load in annotation information
beeAnn <- read.csv("data/Amel_Gene_Names.csv") #Generated from the gff file
antAnn <- read.csv("data/MpharAnn.csv") #This file is from Warner et al 2017 (MBE) and contains other information
colnames(antAnn)[2] = "GeneName"

antGrC = merge(antGrC,antAnn[,c(1,2)],by="Gene",all.x=T)
beeGrC = merge(beeGrC,beeAnn,by="Gene",all.x=T)

allC = rbind(antGrC,beeGrC)
allC$DEcat = factor(allC$DEcat,levels = c("queen","worker","non-DE"))
allC$conBic = as.factor(allC$conBic)
levels(allC$conBic) = c("shared","not shared")

#Make a few big genes have big dots
allC$size=1.9
allC$size[
  grepl("itellogenin",allC$GeneName) | 
  grepl("Smaug",allC$GeneName) | 
  grepl("vasa",allC$GeneName) | 
  grepl("ovo",allC$GeneName)
]=2

#Note -abdomen because we want queen genes to be upregulated in the graph
p1 <- ggplot(allC,aes(x = conn,y=-abdomen))+
  geom_point(aes(fill = DEcat,size=size),pch=21,color="black")+
  scale_fill_manual(values = SexPal[c(1,2,5)],name="DEG type")+
  scale_alpha_continuous(guide="none")+
  ylab("abdominal caste bias")+
  xlab("scaled connectivity in abdominal module")+
  plot2theme+
  facet_grid(. ~ species)+
  theme(panel.spacing = unit(2, "lines"))+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  scale_x_log10()+
  theme(legend.position ="none",
        legend.text = element_text(size=15),
        strip.text = element_text(size=22,face="bold"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        plot.margin = unit(rep(1,4),"cm"),
        strip.background = element_rect(color=NA,fill=NA))+
  annotate("text",x=0.007,y=10,label=c("A","B"),size=16,fontface="bold")

p2 <- ggplot(allC,aes(x=conBic,y=conn))+
  geom_violin(trim=T,position = position_dodge(width=0.8),fill="grey90")+
  geom_jitter(width=0.1,size=1)+
  geom_boxplot(coef=0,size=1,width=0.15,fill="black",color="white",outlier.shape = NA,notch=TRUE,notchwidth = 0.6,position = position_dodge(width=0.8))+
  plot2theme+
  scale_color_manual(values = DE_palette[c(5,9)])+
  ylab("scaled connectivity")+
  facet_grid(. ~ species)+
  theme(panel.spacing = unit(2, "lines"))+
  xlab("presence in abdominal module")+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  theme(legend.position ="none",
        legend.text = element_text(size=15),
        strip.text = element_text(size=22,face="bold"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        plot.margin = unit(c(0,1,1,1),"cm"),
        strip.background = element_rect(color=NA,fill=NA))+
  annotate("text",x=0.62,y=0.9,label=c("C","D"),size=16,fontface="bold")+
  annotate("text",x=1.5,y=0.8,label="***",size=10)
  
  

p <- arrangeGrob(p1,p2,nrow=2)
ggsave(p,file="figures/Fig2.png",width=12,height=12,dpi=300)

###################

#Figure 3

#Identify biased genes and output logFC
extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
}

AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias) #3-5 is head, thorax, abdomen
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

#Label genes by abdominal differential expression conservation
ogg11 = ACUogg
ogg11$abdDE = "non-conserved/non-DE"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "conserved queen" #Note abdomen is the third list
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "conserved worker"

#Combine caste and sex results
FC = merge(AsexRes[[3]][[1]],AcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Mphar",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-conserved/non-DE"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","non-conserved/non-DE"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-conserved/non-DE"]=0.8 #Highlights genes with conserved expression patterns

p1 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal,name = "abdominal caste bias")+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  scale_alpha_continuous(guide="none")+
  main_theme+
  ylab("ant caste bias (queen/worker)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  annotate("text",label="A",size=18,fontface="bold",x = -8,y=8)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",
        legend.text = element_text(size=17),
        legend.title = element_text(size=19,face="bold"))

#Merge caste and sex for honey bees
FC = merge(BsexRes[[3]][[1]],BcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Amel",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-DE/inconsistent"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","non-DE/inconsistent"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-DE/inconsistent"]=0.8

p2 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  annotate("text",label="B",fontface="bold",size=18,x = -8,y=8)+
  ylab("bee caste bias (queen/worker)")+
  xlab("bee sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())

pA <- ggplot()
ggsave(arrangeGrob(p1,p2,nrow=1),file="figures/Fig3ab.png",height=5,width=12)

#compare sex bias between bees and ants
FCb = merge(BsexRes[[3]][[1]],ogg11,by.x="Gene",by.y="gene_Amel")
FC2 = merge(FCb,AsexRes[[3]][[1]],by.x = "gene_Mphar",by.y= "Gene")
FC2$abdDE = factor(FC2$abdDE,levels = c("conserved queen","conserved worker","non-conserved/non-DE"))
FC2$alpha = 0.2
FC2$alpha[FC2$abdDE!="non-conserved/non-DE"]=0.8

p3 <- ggplot(FC2,aes(x = -FC.y,y = -FC.x))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  annotate("text",label="C",fontface="bold",size=18,x = -8,y=8)+
  ylab("bee sex bias (queen/male)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())

#Compare sex bias to drosophila
DmelSC = merge(sexGenes,ogg11,by="gene_Amel")
DmelSC$abdDE = factor(DmelSC$abdDE,levels = c("conserved queen","conserved worker","non-DE/inconsistent"))
DmelSC$alpha = 0.1
DmelSC$alpha[DmelSC$abdDE!="non-DE/inconsistent"]=0.8
DQ = DmelSC[grepl("conserved queen",DmelSC$abdDE),]
DW = DmelSC[grepl("conserved worker",DmelSC$abdDE),]
levels(DmelSC$abdDE) = c("shared queen","shared worker","not caste-biased")

p4 <- ggplot(DmelSC[grepl("shared",DmelSC$abdDE),],aes(x = abdDE,y=-logFC))+ #Female-associated will be positive (intially negative based on DEG direction)
  geom_violin(fill="grey90",trim=FALSE)+
  geom_jitter(width = 0.1,size=0.5,aes(color=abdDE))+
  geom_boxplot(width=0.05,outlier.shape = NA,fill="black",color="black",notch=TRUE,notchwidth = 0.7)+
  plot2theme+
  ylab("fly sex bias (female/male)")+
  ylim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=0,b=0)))+
  xlab("abdominal caste bias")+
  scale_fill_manual(values=SexPal)+
  scale_color_manual(values=SexPal)+
  theme(legend.position="none")+
  annotate("text",label="D",size=18,fontface="bold",x = 0.75,y=6)+
  stat_summary(geom = "crossbar", width=0.035, fatten=0, size=0.7,color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  coord_cartesian(ylim = c(-8,8))+
  theme(legend.position="none")

ggsave(p4,file="figures/Fig3D.png",height=5,width=5,dpi=600)

p <- arrangeGrob(p1,p2,p3,p4,nrow=4)
ggsave(p,file="figures/Fig3.png",height=20,width=5,dpi=300)

#Test for enrichment of genes
binom.test(sum(DmelSC$abdDE=="conserved queen" & DmelSC$logFC < 0,na.rm=T),sum(DmelSC$abdDE=="conserved queen",na.rm=T),0.5)
binom.test(sum(DmelSC$abdDE=="conserved worker" & DmelSC$logFC > 0,na.rm=T),sum(DmelSC$abdDE=="conserved worker",na.rm=T),0.5)

###################

#Figure 4

#Calculate caste/behavior bias as euclidean distance of logFC
euclDist <- function(res){
  cb = apply(res[,-c(1)],1,function(x) sqrt(sum(x^2))/length(x))
  y = res[,colnames(res)!="abdomen"]
  cb_noAbd = apply(y[,-c(1)],1,function(x) sqrt(sum(x^2))/length(x))
  results = data.frame(Gene = res$Gene,cb=cb,cb_noAbd=cb_noAbd)
  return(results)
}

antCB = euclDist(antRes_allstage[[1]]) #This includes larval stages
beeCB = euclDist(beeRes_allstage[[1]])
antSB = euclDist(antSocRes[[1]])
beeSB = euclDist(beeSocRes[[1]])

antCB_expr = euclDist(antRes_allstage[[4]]) #Analogously calculate overall expression
beeCB_expr = euclDist(beeRes_allstage[[4]])
antSB_expr = euclDist(antSocRes[[4]])
beeSB_expr = euclDist(beeSocRes[[4]])

#Get results together
antCB = merge(antCB,antCB_expr,by="Gene")
beeCB = merge(beeCB,beeCB_expr,by="Gene")
antSB = merge(antSB,antSB_expr,by="Gene")
beeSB = merge(beeSB,beeSB_expr,by="Gene")

colnames(antCB) = colnames(beeCB) = colnames(antSB) = colnames(beeSB) =c("Gene","cb","cb_noAbd","expr","expr_noAbd")

#Expression is negatively correlated to bias
cor.test(antCB$cb,antCB$expr,method="spearman")
cor.test(antSB$cb,antSB$expr,method="spearman")
cor.test(beeCB$cb,beeCB$expr,method="spearman")
cor.test(beeSB$cb,beeSB$expr,method="spearman")

antCB$type=beeCB$type="caste"
antSB$type=beeSB$type="behavior"
aC = rbind(antCB,antSB)
bC = rbind(beeCB,beeSB)

aC = merge(aC,Aps,by="Gene")
bC = merge(bC,Bps,by="Gene")
aC = merge(antConn,aC,by="Gene")
bC = merge(beeConn,bC,by="Gene")

#Scale connectivity to top value
aC$kTotal = aC$kTotal/max(aC$kTotal)
bC$kTotal = bC$kTotal/max(bC$kTotal)

#Add dn/ds data
antD <- read.table("results/mphar_sinv_dnds.txt",head=T)
beeD <- read.table("results/amel_cerana_dnds.txt",head=T)

aC = merge(aC,antD,by="Gene",all.x=T)
bC = merge(bC,beeD,by="Gene",all.x=T)


levels(aC$psName)[1] = levels(bC$psName)[1]= "ancient"

#Make ordinal variable for statistical analyses
aC$psOrd = ordered(aC$psName)
bC$psOrd = ordered(bC$psName)

aS = aC[aC$type=="behavior",]
bS = bC[bC$type=="behavior",]
aC = aC[aC$type=="caste",]
bC = bC[bC$type=="caste",]

cor.test(aC$cb,aC$dN_dS,method="spearman")
cor.test(bC$cb,bC$dN_dS,method="spearman")
cor.test(aS$cb,aS$dN_dS,method="spearman")
cor.test(bS$cb,bS$dN_dS,method="spearman")

cor.test(aC$cb,aC$kTotal,method="spearman")
cor.test(bC$cb,bC$kTotal,method="spearman")
cor.test(aS$cb,aS$kTotal,method="spearman")
cor.test(bS$cb,bS$kTotal,method="spearman")

#Make plots for relationship between caste bias and phylostrata, connectivity, and dn/ds
p1 <- ggplot(aC,aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("ant")+
  xlab("estimated evolutionary age")+
  coord_cartesian(ylim=c(0,2.5))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  annotate("text",label="A",size=18,fontface="bold",x = 5.5,y=2.2)+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")


p2 <- ggplot(bC[!is.na(bC$psName),],aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("honey bee")+
  xlab("estimated evolutionary age")+
  coord_cartesian(ylim=c(0,2.5))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  annotate("text",label="B",size=18,fontface="bold",x = 5.5,y=2.2)+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")

p3 <- ggplot(aC,aes(x = kTotal,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("scaled network connectivity")+
  ggtitle("ant")+
  ylab("overall caste bias")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  annotate("text",label="C",size=18,fontface="bold",x = 0.4,y=2.2)+
  theme(legend.position="none")

p4 <- ggplot(bC,aes(x = kTotal,y=cb))+
  geom_hex(bins=70)+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  scale_fill_gradient(low = "blue",high="red")+
  xlab("scaled network connectivity")+
  ylab("overall caste bias")+
  ggtitle("honey bee")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  annotate("text",label="D",size=18,fontface="bold",x = 0.4,y=2.2)+
  theme(legend.position="none")

p5 <- ggplot(aC,aes(x = dN_dS,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,2.5)+
  xlim(0,1.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("dN/dS")+
  ggtitle("ant")+
  ylab("overall caste bias")+
  annotate("text",label="E",size=18,fontface="bold",x = 1.28,y=2.2)+
  #scale_x_log10(breaks = c(0.01,0.1,1))+
  theme(legend.position="none")

p6 <- ggplot(bC,aes(x = dN_dS,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  xlim(0,1.5)+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("dN/dS")+
  ggtitle("honey bee")+
  ylab("overall caste bias")+
  #scale_x_log10(breaks = c(0.01,0.1,1))+
  annotate("text",label="F",size=18,fontface="bold",x = 1.28,y=2.2)+
  theme(legend.position="none")

p <- arrangeGrob(p1,p2,p3,p4,p5,p6,nrow=3)
ggsave(p,file = "figures/Fig4.png",height=15,width=10,dpi=300)

ggsave(arrangeGrob(p5,p6,nrow=1),file="~/Downloads/Fig4ef.png",height=5,width=10,dpi=300)
ggsave(arrangeGrob(p1,p2,nrow=1),file="~/Downloads/Fig4ab.png",height=5,width=10,dpi=300)

ggsave(arrangeGrob(p3,p4,nrow=1),file="~/Downloads/Fig4cd.png",height=5,width=10,dpi=300)


########
#Figure S8 - behavior bias vs evolutionary statistics

p1 <- ggplot(aS,aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("ant")+
  xlab("estimated evolutionary age")+
  coord_cartesian(ylim=c(0,4))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall behavior bias")+
  annotate("text",label="A",size=18,fontface="bold",x = 5.5,y=3.5)+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")


p2 <- ggplot(bS[!is.na(bS$psName),],aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("honey bee")+
  xlab("estimated evolutionary age")+
  coord_cartesian(ylim=c(0,4))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall behavior bias")+
  annotate("text",label="B",size=18,fontface="bold",x = 5.5,y=3.5)+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")

p3 <- ggplot(aS,aes(x = kTotal,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,4)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("scaled network connectivity")+
  ggtitle("ant")+
  ylab("overall behavior bias")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  annotate("text",label="C",size=18,fontface="bold",x = 0.4,y=3.5)+
  theme(legend.position="none")

p4 <- ggplot(bS,aes(x = kTotal,y=cb))+
  geom_hex(bins=70)+
  plot2theme+
  ylim(0,4)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  scale_fill_gradient(low = "blue",high="red")+
  xlab("scaled network connectivity")+
  ylab("overall behavior bias")+
  ggtitle("honey bee")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  annotate("text",label="D",size=18,fontface="bold",x = 0.4,y=3.5)+
  theme(legend.position="none")

p5 <- ggplot(aS,aes(x = dN_dS,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,4)+
  xlim(0,1.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("dN/dS")+
  ggtitle("ant")+
  ylab("overall behavior bias")+
  annotate("text",label="E",size=18,fontface="bold",x = 1.28,y=3.5)+
  theme(legend.position="none")

p6 <- ggplot(bS,aes(x = dN_dS,y=cb))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  xlim(0,1.5)+
  ylim(0,4)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("dN/dS")+
  ggtitle("honey bee")+
  ylab("overall behavior bias")+
  annotate("text",label="F",size=18,fontface="bold",x = 1.28,y=3.5)+
  theme(legend.position="none")

p <- arrangeGrob(p1,p2,p3,p4,p5,p6,nrow=3)
ggsave(p,file = "figures/FigS8.png",height=15,width=10,dpi=300)

############

#Make table of parital correlations (Table S11)

#Add tau for honey bees
tau <- read.csv("results/bee_tau.csv")
bC2 = merge(tau,bC,by="Gene")
bS2 = merge(tau,bS,by="Gene")

pCB <- function(data,tau=FALSE){
  data$psNum = as.numeric(data$psName)
  data = data[!is.na(data$dN_dS),]
  tests <- list(data$kTotal,data$dN_dS,data$psNum)
  var = c("connectivity","dN_dS","evolutionary age")
  if (tau){
    tests <- list(data$kTotal,data$dN_dS,data$psNum,data$tau)
    var = c("connectivity","dN_dS","evolutionary age","tau")
  }
  cb <- list(data$cb,data$cb_noAbd)
  expr <- list(data$expr,data$expr_noAbd)
  aK = c("yes","no")
  results <- ldply(lapply(c(1:length(var)),function(i){
    ldply(lapply(c(1,2),function(j){
      data.frame(Abd = aK[j],variable=var[i],pcor.test(cb[[j]],tests[[i]],expr[[j]],method="spearman"))
    }))
  }))
  return(results)
}

a1 <- pCB(aC)
b1 <- pCB(bC2,tau=T)
a2 <- pCB(aS)
b2 <- pCB(bS2,tau=T)

a1$species = a2$species="ant"
b1$species = b2$species = "honey bee"
a1$type=b1$type="caste"
a2$type=b2$type="behavior"

pCor <- do.call(rbind,list(a1,b1,a2,b2))
pCor = pCor[,c(9,10,1:4)]
pCor$variable = as.character(pCor$variable)
pCor$variable[pCor$variable=="dN_dS"]="dN/dS"
colnames(pCor) = c("species","comparison","abdomen included?","variable tested","Spearman rho","P-value")
pCor$`Spearman rho` = round(pCor$`Spearman rho`,3)
pCor$`P-value` = signif(pCor$`P-value`,3)

tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(pCor, rows=NULL, theme=tt)
ggsave(tbl,file = "figures/Table_pcor.png",height=12,width=10,dpi=300)
################

#Figure S1
#Correlation of log fold change across development
lfcCor <- function(antD,beeD){
  nStage = ncol(antD) - 1 
  antD <- merge(antD,ACUogg,by.x = "Gene",by.y = "gene_Mphar")
  beeD <- merge(beeD,ACUogg,by.x = "Gene",by.y = "gene_Amel")
  antD = antD[antD$OGG %in% beeD$OGG,]
  beeD = beeD[beeD$OGG %in% antD$OGG,]
  antD = antD[order(antD$OGG),]
  beeD = beeD[order(beeD$OGG),]
  d = data.frame(stage = colnames(antD)[c(2:(nStage+1))])
  dAbs = data.frame(stage = colnames(antD)[c(2:(nStage+1))])
  for (i in 1:nStage){
    t = cor.test(antD[,i+1],beeD[,i+1])
    d[i,2] = t$estimate
    d[i,3] = t$conf.int[1]
    d[i,4] = t$conf.int[2]
    d[i,5] = t$p.value
  }
  colnames(d) = c("Stage","cor","c1","c2","pval")
  return(d)
}

CasteCor <- lfcCor(antRes[[1]],beeRes[[1]])
CasteCor_allStage <- lfcCor(antRes_allstage[[1]],beeRes_allstage[[1]])
BehavCor <- lfcCor(antSocRes[[1]],beeSocRes[[1]])

CasteCor$Stage = as.character(CasteCor[[1]]$Stage)
CasteCor$Stage[1] = "larva*"
CasteCor_allStage$Stage = c("L2","L3","L4","L5","pupa","head","thorax","abdomen")
BehavCor$Stage = c("head","thorax","abdomen")
CasteCor$Stage = factor(CasteCor$Stage,levels = CasteCor$Stage)
CasteCor_allStage$Stage = factor(CasteCor_allStage$Stage,levels = CasteCor_allStage$Stage)
BehavCor$Stage = factor(BehavCor$Stage,levels = BehavCor$Stage)

#Put caste figures together
pl <- lapply(list(CasteCor,CasteCor_allStage,BehavCor), function(x){
  ggplot(x,aes(x = Stage,y=cor))+
    geom_bar(stat="identity")+main_theme+
    xlab("stage/tissue")+
    geom_errorbar(aes(ymin=c1,ymax=c2),width=0.2)+
    ylab("inter-species correlation of\nqueen/worker log2 fold-change")+
    theme(plot.margin=unit(c(0.5,1.5,0.5,0.5),"cm"),
          axis.text.x = element_text(hjust=0,angle=-25))
})

pl[[3]] = pl[[3]]+xlab("tissue")+
  ylab("correlation of nurse/forager log2 fold-change")

p <- arrangeGrob(pl[[1]]+annotate("text",x=1,y=0.25,label="A",size=18,fontface="bold"),
                 pl[[2]]+annotate("text",x=1.3,y=0.25,label="B",size=18,fontface="bold"),
                 ncol=2)
ggsave(p,file="figures/FigS1.png",height=8,width=13.5,dpi=300)

ggsave(pl[[3]],file="figures/FigS2.png",height=8,width=7,dpi=300)

############


############

#Figure S3

#Merge phylostrata information with DEG information
Aps$species = "M. pharaonis"
Bps$species = "A. mellifera"
Aps_FC = merge(Aps,antRes[[1]],by="Gene")
Bps_FC = merge(Bps,beeRes[[1]],by="Gene")
Bps_FC = Bps_FC[!is.na(Bps_FC$psName),]

#Adjust logFC values to medians
for (i in 7:11){
  Aps_FC[,i] = Aps_FC[,i] - median(Aps_FC[,i])
  Bps_FC[,i] = Bps_FC[,i] - median(Bps_FC[,i])
}

levels(Aps_FC$psName)[1:4] = levels(Bps_FC$psName)[1:4] = c("ancient","insect","hymenopteran","aculeate")

pl <- lapply(list(Aps_FC,Bps_FC),function(x){
  pm = melt(x,id.vars = c("Gene","ODBgene","OGGacu","ps","psName","species"))
  ggplot(pm,aes(x = variable,y=-value,fill=psName))+
    geom_boxplot(notch=TRUE,outlier.shape = NA)+
    coord_cartesian(ylim = c(-6,6))+
    main_theme+
    scale_fill_manual(values = rev(mypalette),name = "phylostrata")+
    xlab("stage/tissue")+
    ylab("log2 fold-change\n(queen/worker)")+
    geom_hline(yintercept =0,linetype="dashed")+
    theme(legend.position="right",
          #legend.key.size = unit(0.6,"cm"),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15,face="bold"),
          legend.background = element_blank(),
          plot.margin = unit(rep(1,4),"cm"))
})

p <- arrangeGrob(pl[[1]]+annotate("text",x=1.4,y=-5,label="A. ant",size=11,fontface="bold"),
                 pl[[2]]+annotate("text",x=1.6,y=-5,label="B. honey bee",size=11,fontface="bold"),
                 ncol=1)

ggsave(p,file="~/GitHub/devnetwork/figures/FigS2.png",height=12,width=10,dpi=300)

############

#Figure S4
antD = rownames(antDevel2)[antDevel2$FDR < 0.1]
beeD = rownames(beeDevel2)[beeDevel2$FDR < 0.1]

ACUogg$beeD = ACUogg$antD = 0
ACUogg$beeD[ACUogg$gene_Amel %in% beeD] = 1
ACUogg$antD[ACUogg$gene_Mphar %in% antD] = 1

a = as.data.frame(rbind(c(length(antD),sum(ACUogg$antD),sum(ACUogg$antD & ACUogg$beeD)),
                        c(length(beeD),sum(ACUogg$beeD),sum(ACUogg$antD & ACUogg$beeD))))

a$species = c("ant","honey bee")
colnames(a)[1:3] = c("develTotal","devel_withOGG","devel_bothSpecies")
a$devTRG = a$develTotal - a$devel_withOGG
a$OGGns = a$devel_withOGG - a$devel_bothSpecies
aM = melt(a[,c(3:6)],id.vars = "species")
aM$variable = factor(aM$variable,levels = levels(aM$variable)[c(2,3,1)])
levels(aM$variable) = c(" no ortholog   "," not shared developmental   "," shared developmental   ")

p1 <- ggplot(aM,
              aes(x = species, y = value, fill = variable))+
  geom_bar(stat="identity",color="black")+
  plot2theme+
  xlab("species")+
  scale_fill_manual(values = DE_palette[c(3,7,9)])+
  ylab("number of\ndevelopmental DEGs")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(panel.spacing = unit(2, "lines"))

ggsave(p1,file = "figures/FigSdevel.png",height=8,width=8,dpi=300)
#############

#Figure S4

#Calculate correlation of caste and sex bias
FCcor <- function(test1,test2){
  FC = merge(test1[[1]],test2[[1]],by = "Gene")
  r = cor.test(FC$FC.x,FC$FC.y)
  return(c(r$estimate,c1=r$conf.int[1],c2=r$conf.int[2]))
}

#Calculate correlation of caste/sex bias across tissues
Acor <- ldply(lapply(seq(1,3),function(i){
  FCcor(AcasteRes[[i]],AsexRes[[i]])
}))

Bcor <- ldply(lapply(seq(1,3),function(i){
  FCcor(BcasteRes[[i]],BsexRes[[i]])
}))

Acor$tissue = Bcor$tissue = c("head","thorax","abdomen")
Acor$species = "ant"
Bcor$species = "honey bee"
AllCor = rbind(Acor,Bcor)
AllCor$tissue=factor(AllCor$tissue,levels = c("head","thorax","abdomen"))

p <- ggplot(AllCor,aes(x=species,y=cor,fill=tissue))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(.9),width=0.2)+
  main_theme+
  ylim(0,1)+
  scale_fill_manual(values=c("#edf8b1","#7fcdbb","#2c7fb8"))+
  ylab("correlation of caste and sex bias")+
  theme(legend.title = element_text(size=18,face="bold"),
        legend.position = c(0.15,0.8))
ggsave(p,file = "figures/FigS4.png",height=6,width=8)

#Figure S5
#Summarize number of times DE for caste-associated genes
sumDE <- function(dfDE,type1,type2){
  dfDE$numQueen = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == type1))
  dfDE$numWorker = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == type2))
  d = table(dfDE$numQueen,dfDE$numWorker)
  m = melt(d)
  colnames(m)[c(1,2)] = c(type1,type2)
  return(m)
}


m1 = sumDE(antRes[[2]],"queen","worker")
m2 = sumDE(beeRes[[2]],"queen","worker")

#Add these since there are no genes DE all five times in Apis 
m2E = t(sapply(seq(0,5),function(i) c(queen=i,worker=5,value=0)))
m2Eb = t(sapply(seq(0,5),function(i) c(queen=5,worker=i,value=0)))
m2E = rbind(m2Eb,m2E)
m2 = rbind(m2,m2E)

m1$species = "ant"
m2$species = "honey bee"
mA = rbind(m1,m2)


#Create heatmap of differential expression (number of times DE for queens and workers)
p <- ggplot(mA,aes(x=queen,y=worker))+
  geom_tile(aes(fill = value))+
  facet_grid(. ~ species)+
  scale_fill_gradient(name = "number of genes",trans = "log",
                      breaks = c(1,10,100,1000),
                      limits = c(1,10000),
                      labels = c(1,10,100,1000))+
  geom_text(aes(x = queen,y = worker,label = value),color="white")+
  main_theme+
  scale_y_continuous(name = "number of times worker-biased",
                     breaks  = seq(0,5),
                     expand = c(0,0))+
  scale_x_continuous(name = "number of times queen-biased",
                     breaks = seq(0,5),
                     expand = c(0,0))+
  theme(legend.position = "right",
        axis.line=element_line(color="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        panel.border = element_rect(size = 1, color = "black",fill = NA)) 
ggsave(p,file = "figures/FigS5.png",height=8,width=10)

###############

#Figure S6 - compare caste bias between species
cbBias = merge(beeCB,ACUogg,by.x="Gene",by.y="gene_Amel") #ACUogg has ortholog definitions
cbBias = merge(cbBias,antCB,by.x="gene_Mphar",by.y="Gene")
cbBias = cbBias[!is.na(cbBias$OGGacu),] #Remove genes without 1-1 orthologs
socBias = merge(beeSB,ACUogg,by.x="Gene",by.y="gene_Amel")
socBias = merge(socBias,antSB,by.x="gene_Mphar",by.y="Gene")
socBias = socBias[!is.na(socBias$OGGacu),]
p3 <- ggplot(cbBias,aes(x=cb.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  annotate("text",size=16,fontface="bold",x=2,y=1.75,label="A")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  main_theme+xlab("honey bee caste bias")+ylab("ant caste bias")

p4 <- ggplot(socBias,aes(x=cb.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",size=16,fontface="bold",x=4,y=3.5,label="B")+
  main_theme+xlab("honey bee behavioral bias")+ylab("ant behavior bias")

p <- arrangeGrob(p3,p4,nrow=1)
ggsave(p,file="figures/FigS6.png",width=12,height=7)
ggsave(p,file="figures/FigS6.png",width=12,height=6)


###############

#Figure S7 compare caste and behavior bias within species
antBias = merge(antCB,antSB,by="Gene")
beeBias = merge(beeCB,beeSB,by="Gene")
p3 <- ggplot(antBias,aes(x=cb.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  xlim(0,2.5)+
  ylim(0,4)+
  annotate("text",size=16,fontface="bold",x=2,y=3.5,label="A")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  main_theme+xlab("ant caste bias")+ylab("ant behavior bias")

p4 <- ggplot(beeBias,aes(x=cb.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  xlim(0,2.5)+
  ylim(0,4)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",size=16,fontface="bold",x=2,y=3.5,label="B")+
  main_theme+xlab("honey bee caste bias")+ylab("honey bee behavior bias")

p <- arrangeGrob(p3,p4,nrow=1)
ggsave(p,file="figures/FigS7.png",width=12,height=7)


#Figure S9 - compare bias to tissue specificity (tau) in honey bees
tau <- read.csv("results/bee_tau.csv")
bT = merge(tau,beeCB,by="Gene")

p1 <- ggplot(bT,aes(x=tau,y=cb))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(se=FALSE,method="lm",color="red")+
  ylim(0,2.5)+
  annotate("text",x = 0.3,y = 2,label="A",size=16,fontface="bold")+
  main_theme+xlab("tissue specificity")+ylab("honey bee caste bias")

bT = merge(tau,beeSB,by="Gene")

p2 <- ggplot(bT,aes(x=tau,y=cb))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(se=FALSE,method="lm",color="red")+
  ylim(0,5)+
  annotate("text",x = 0.3,y = 4,label="B",size=16,fontface="bold")+
  main_theme+xlab("tissue specificity")+ylab("honey bee behavior bias")
p <- arrangeGrob(p1,p2,nrow=1)
ggsave(p,file="figures/FigS9.png",height=7,width=12)

##############


###############

#Fig S10
load("results/PlaidResults.RData")
freq_set <- function(data,maxMod,factor,tpm,type){
  res <- lapply(data,function(x){
    checkMod = seq(1:nrow(x@NumberxCol))[rowSums(x@NumberxCol) <= maxMod]
    r = list()
    modI = 0
    for (mod in checkMod){
      qp = x@NumberxCol
      samp = factor$sample[x@NumberxCol[mod,]]
      if (sum(grepl(type,samp))==3){
        modI = modI+1
        genes = rownames(tpm)[x@RowxNumber[,mod]]
        r[[modI]]=list(samples=samp,genes=genes)
      }
    }
    if (length(r) > 1){
      rnew = r[[1]]
      rnew$genes = unique(c(r[[2]]$genes,r[[1]]$genes)) #combine genes from queen abdomen module
      r = rnew
    }
    return(r)
  })
  return(res)
}

#Calculate number of significant biclusters for a given set
num_set <- function(set){
  nums = lapply(names(set),function(x){
    length(set[[x]][lapply(set[[x]],length) > 0])
  })
  numDF = data.frame(Type = names(set),Freq=unlist(nums))
}

commonGenes <- function(QGres,tpm){
  df = data.frame(Gene=rownames(tpm),KeepNum=0)
  num = 0
  for (res in QGres){
    if (length(res) == 1){
      num = num+1
      df$KeepNum[df$Gene %in% res[[1]]$genes] = df$KeepNum[df$Gene %in% res[[1]]$genes] + 1
    }
  }
  df$KeepNum = df$KeepNum/num
  return(df)
}

samp_types = sapply(factorA$sample,function(x) gsub(".*_","",x) %>% gsub("\\.[0-9]","",.)) %>%
  unique()

#Output significant biclusters for each sample type
type_bicAmel <- lapply(samp_types,function(x) freq_set(beePl,6,factorB,beeT,x))
type_bicMphar <- lapply(samp_types,function(x) freq_set(antPl,6,factorA,antT,x))
names(type_bicAmel) = names(type_bicMphar) = samp_types
typeNum_Amel = num_set(type_bicAmel)
typeNum_Mphar = num_set(type_bicMphar)

antQG <- type_bicMphar[["AQG"]]
beeQG <- type_bicAmel[["AQG"]]
antQG = antQG[lapply(antQG,length) > 0]
beeQG = beeQG[lapply(beeQG,length) > 0]

aG <- commonGenes(antQG,antT)
bG <- commonGenes(beeQG,beeT)

p1 <- ggplot(aG,aes(x=KeepNum))+
  geom_histogram(bins = 20,fill="white",color="black")+
  plot2theme+
  ylim(0,5000)+
  xlab("frequency of presence in bicluster")+
  ylab("number of genes")+
  annotate("text",size=12,fontface="bold",x=0.8,y=4500,label="A. ant")
p2 <- ggplot(bG,aes(x=KeepNum))+
  geom_histogram(bins = 20,fill="white",color="black")+
  plot2theme+
  ylim(0,9000)+
  xlab("frequency of presence in bicluster")+
  ylab("number of genes")+
  annotate("text",size=12,fontface="bold",x=0.65,y=4500*1.8,label="B. honey bee")

p <- arrangeGrob(p1,p2,nrow=1)
ggsave(p,file="figures/S3.png",height=6,width=12,dpi=300)

ant$Gene = rownames(ant)
ant2 = merge(ant,ACUogg,by.x="Gene",by.y="gene_Mphar")
bee$Gene = rownames(bee)
bee2 = merge(ACUogg,bee,by.y="Gene",by.x="gene_Amel")
all = merge(ant2,bee2,by="OGGacu")
all = all[,-c(1,2,94)]
all = all[,c(91,92,1:90,93:179)]
colnames(all)[1] = "gene_Amel"
factorA$species = "ant"
factorB$species = "bee"
factors = rbind(factorA,factorB)
write.csv(factors,file="data/allFactor.csv",row.names=F)
write.csv(all,file="data/OrthoExpr.csv",row.names=F)
