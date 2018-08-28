load("~/GitHub/devnetwork/results/DEtests.RData")
load("~/GitHub/devnetwork/results/collectedPhylo.RData")
setwd("~/GitHub/devnetwork/src/analysis_R/")
load("../../results/PlaidResults.RData")

TGmap <- read.table("~/GitHub/devnetwork/phylo_results/TGmap_Amel.txt")
TNmap <- as.data.frame(fread("~/GitHub/devnetwork/data/AmelTranName.txt",sep="~",header=FALSE))

AmelName <- merge(TGmap,TNmap,by.x = "V2",by.y = "V1")[,c(2,3)]
colnames(AmelName) = c("Gene","GeneName")
AmelName$GeneName = gsub(" isoform X[0-9]","",AmelName$GeneName)
aName = AmelName[!duplicated(AmelName$Gene),]

library(cowplot)
library(RColorBrewer)
library(VennDiagram)
library(plyr)
library(data.table)
library(zoo)
library(reshape2)
library(EnvStats)
library(grid)
library(gridExtra)
DE_palette = brewer.pal(9,"Blues")
mypalette <- brewer.pal(6,"OrRd")
mypalette2 <- brewer.pal(6,"Blues")
SexPal = c("firebrick2","slateblue4","yellow2","darkgoldenrod2","gray94")
specPal = c("gold1","#b33c00")

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right","top")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none",
                                            axis.line.x = element_line(color='black'),
                                            axis.line.y = element_line(color='black'),
                                            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")))
  gl <- c(gl, ncol = ncol, nrow = 1)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)),
                     "top" = arrangeGrob(legend,
                                         do.call(arrangeGrob, gl),
                                         ncol = 1,
                                         heights = unit.c(lheight,unit(1, "npc") - lheight))
                     
  )
  
  return(combined)
  
}

grid_arrange_shared_legend2 <- function(..., ncol = length(list(...)), nrow = nrow, position = c("bottom", "right","top")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none",
                                            axis.line.x = element_line(color='black'),
                                            axis.line.y = element_line(color='black'),
                                            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)),
                     "top" = arrangeGrob(legend,
                                         do.call(arrangeGrob, gl),
                                         ncol = 1,
                                         heights = unit.c(lheight,unit(1, "npc") - lheight))
                     
  )
  
  return(combined)
  
}
grid_arrange_shared_legend <- function(...,ncol) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none")),list(ncol=ncol)),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
euclDist <- function(res){
  cb = apply(res[,-c(1)],1,function(x) sqrt(sum(x^2))/length(x))
  cb_noAbd = apply(res[,-c(1,ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_noAdult = apply(res[,-c(1,(ncol(res) - 2):ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_larva = apply(res[,-c(1,(ncol(res) - 3):ncol(res))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_adult = apply(res[,-c(1:(ncol(res) - 3))],1,function(x) sqrt(sum(x^2))/length(x))
  cb_abd = apply(as.data.frame(res[,-c(1:(ncol(res) - 1))]),1,function(x) sqrt(sum(x^2))/length(x))
  results = data.frame(Gene = res$Gene,cb=cb,cb_noAbd=cb_noAbd,cb_noAdult=cb_noAdult,cb_larva=cb_larva,cb_adult = cb_adult,cb_abd=cb_abd)
  return(results)
}
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


#Get coefficient of variation for a particular sample type
getCV <- function(lev,factor,tpm){
  if (sum(factor$meta==lev)<2){
    return(rep(0,nrow(tpm)))
  }
  samples = factor$sample[factor$meta==lev]
  df = tpm[,samples]
  cv = apply(df,1,function(x) sd(x)/mean(x))
  return(cv)
}

#Average coefficient of variation across all sample types
averageCV <- function(factor,tpm){
  allCv = lapply(levels(factor$meta),function(x) getCV(x,factor,tpm))
  a = as.data.frame(do.call(cbind,allCv))
  names(a) = levels(factor$meta)
  meanCv = apply(a,1,function(x) mean(x,na.rm=TRUE))
  return(list(meanCv,a))
}


geomMean <- function(res){
  cb = apply(res[,-c(1)],1,function(x) gm_mean(x))
  cb_noAbd = apply(res[,-c(1,ncol(res))],1,function(x) gm_mean(x))
  cb_noAdult = apply(res[,-c(1,(ncol(res) - 2):ncol(res))],1,function(x) gm_mean(x))
  cb_larva = apply(res[,-c(1,(ncol(res) - 3):ncol(res))],1,function(x) gm_mean(x))
  cb_adult = apply(res[,-c(1:(ncol(res) - 3))],1,function(x) gm_mean(x))
  cb_abd = apply(as.data.frame(res[,-c(1:(ncol(res) - 1))]),1,function(x) gm_mean(x))
  results = data.frame(Gene = res$Gene,cb=cb,cb_noAbd=cb_noAbd,cb_noAdult=cb_noAdult,cb_larva=cb_larva,cb_adult = cb_adult,cb_abd=cb_abd)
  return(results)
}

modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}



#generate factors df based on counts fie
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

#identify a QG-specific module
idQG <- function(data,maxMod,factor,tpm){
  res <- lapply(data,function(x){
    checkMod = seq(1:nrow(x@NumberxCol))[rowSums(x@NumberxCol) <= maxMod]
    r = list()
    modI = 0
    for (mod in checkMod){
      qp = x@NumberxCol
      samp = factor$sample[x@NumberxCol[mod,]]
      if (sum(grepl("QG",samp))==3){
        modI = modI+1
        genes = rownames(tpm)[x@RowxNumber[,mod]]
        r[[modI]]=list(samples=samp,genes=genes)
      }
    }
    return(r)
  })
  return(res)
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

#Caste-associated genes are also sex-biased
extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
}


#filter out based on counts per million reads
filterLowly <- function(counts,cut){ 
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= (ncol(counts)/2)
  return(counts[keep,])
}


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


beeT <- read.table("../../data/bees.tpm.txt",header=TRUE)
antT <- read.table("../../data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

bee <- read.table("~/GitHub/devnetwork/data/bees.counts_edit.txt",header=TRUE)
ant <- read.table("~/GitHub/devnetwork/data/ants.counts_edit.txt",header=TRUE)
bee = modifyDF(bee)
ant = modifyDF(ant)

#Filter out lowly expressed genes
ant <- filterLowly(ant,1)
bee <- filterLowly(bee,1)

#Make factors
factorA <- genFactor(ant)
factorB <- genFactor(bee)

############
###Fig 1a
############
#Compare DE definition at each stage
antRes[[2]]$ortholog_found = antRes[[2]]$OGG_found = FALSE
antRes[[2]]$ortholog_found[antRes[[2]]$Gene %in% AllPS$Gene.x] = TRUE
antRes[[2]]$OGG_found[antRes[[2]]$Gene %in% ACUogg$gene_Mphar] = TRUE
aM = melt(antRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
beeRes[[2]]$ortholog_found = beeRes[[2]]$OGG_found = FALSE
beeRes[[2]]$ortholog_found[beeRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeRes[[2]]$OGG_found[beeRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
aMcaste = aM #saving for later
bMcaste = bM
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

aD$DEboth = bD$DEboth = allD$DEboth
aD$OGG = aD$OGG - aD$DEboth
bD$OGG = bD$OGG - bD$DEboth

aDM = melt(aD,id.vars = "variable")
bDM = melt(bD,id.vars = "variable")
colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")
aDM$species = "ant"
bDM$species = "bee"

d = rbind(aDM,bDM)
d$species=as.factor(d$species)
levels(d$species) = c("M. pharaonis","A. mellifera")
levels(d$DEtype) = c("NDE","no ortholog","paralogs present","ortholog present,\nnon-conserved caste-bias","ortholog present,\nconserved caste-bias")
levels(d$stage)[1] = "larva*"

p1 <- ggplot(d[d$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  ylim(0,5500)+
  facet_grid(. ~ species)+
  xlab("stage/tissue")+
  scale_fill_manual(values = DE_palette[c(3,5,7,8)])+
  ylab("number of caste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.margin = margin(0,2,2,2,"cm"))

levels(d$species) = c("ant","bee")
p1m <- ggplot(d[d$DEtype!="NDE",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  plot2theme+
  ylim(0,5500)+
  facet_grid(. ~ species)+
  xlab("stage/tissue")+
  scale_fill_manual(values = DE_palette[c(3,5,7,9)])+
  ylab("number of\ncaste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=22,face="bold"),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA,fill=NA),
        plot.margin = margin(0.5,2,0.5,0.5,"cm"))+
  theme(panel.spacing = unit(2, "lines"))

p2 <- ggplot(d[d$DEtype!="NDE" & d$species=="ant",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  scale_fill_manual(values = DE_palette[c(3,5,7,9)])+
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


p3 <- ggplot(d[d$DEtype!="NDE" & d$species=="bee",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  scale_fill_manual(values = DE_palette[c(3,5,7,9)])+
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

pQ <- ggdraw()+
  draw_plot(p1m+
              theme(legend.text = element_text(size=15),
                    legend.key.width = unit(1,"cm"),
                    legend.position = "bottom"))+
  draw_plot(p2,x=0.15,y=0.69,height=0.18,width=0.18)+
  draw_plot(p3,x=0.57,y=0.69,height=0.18,width=0.18)

ggsave(pQ,file = "~/GitHub/devnetwork/figures/fig1a.png",height=6,width=11,dpi=300)

############
###Fig 1b
############
#Compare DE definition at each stage
antSocRes[[2]]$ortholog_found = antSocRes[[2]]$OGG_found = FALSE
antSocRes[[2]]$ortholog_found[antSocRes[[2]]$Gene %in% AllPS$Gene.x] = TRUE
antSocRes[[2]]$OGG_found[antSocRes[[2]]$Gene %in% ACUogg$gene_Mphar] = TRUE
aM = melt(antSocRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
aMsoc = aM
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="non-DE"),
           no_ortholog = sum(value!="non-DE" & !ortholog_found),
           dup = sum(value!="non-DE" & ortholog_found & !OGG_found),
           OGG = sum(value!="non-DE" & OGG_found))
beeSocRes[[2]]$ortholog_found = beeSocRes[[2]]$OGG_found = FALSE
beeSocRes[[2]]$ortholog_found[beeSocRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeSocRes[[2]]$OGG_found[beeSocRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeSocRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bMsoc = bM
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="non-DE"),
           no_ortholog = sum(value!="non-DE" & !ortholog_found),
           dup = sum(value!="non-DE" & ortholog_found & !OGG_found),
           OGG = sum(value!="non-DE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="non-DE" & value != "non-DE"))

aD$DEboth = bD$DEboth = allD$DEboth
aD$OGG = aD$OGG - allD$DEboth
bD$OGG = bD$OGG - allD$DEboth

aDM = melt(aD,id.vars = "variable")
bDM = melt(bD,id.vars = "variable")
colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")
aDM$species = "ant"
bDM$species = "bee"

d = rbind(aDM,bDM)
d$species=as.factor(d$species)
levels(d$species) = c("M. pharaonis","A. mellifera")
levels(d$DEtype) = c("NDE","no ortholog","paralogs present","ortholog present,\nnon-conserved role-bias","ortholog present,\nconserved role-bias")
levels(d$stage) = c("head","thorax","abdomen")

p1 <- ggplot(d[d$DEtype!="NDE",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  xlab("tissue")+
  ylim(0,4500)+
  facet_grid(. ~ species)+
  scale_fill_manual(values = DE_palette)+
  ylab("number of behavior-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.margin = margin(0,2,2,2,"cm"))

p2 <- ggplot(d[d$DEtype!="NDE" & d$species=="M. pharaonis",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  scale_fill_manual(values = DE_palette)+
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


p3 <- ggplot(d[d$DEtype!="NDE" & d$species=="A. mellifera",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  scale_fill_manual(values = DE_palette)+
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
  draw_plot(p1+
              theme(legend.text = element_text(size=13),
                    legend.key.width = unit(1,"cm")))+
  draw_plot(p2,x=0.2,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.58,height=0.18,width=0.18)


psO = merge(AllPS_sum,ACUogg,by="OGGacu")
psO$a2 = "non-DE/non-conserved"
#psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen!="non-DE"] | 
 #        psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen!="non-DE"]] = "non-conserved"
#psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen!="non-DE"] & 
#         psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen!="non-DE"]] = "flipped"

psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="worker"] & 
         psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen=="worker"]] = "conserved worker"
psO$a2[psO$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="queen"] & 
         psO$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen=="queen"]] = "conserved queen"

levels(psO$psName)[1] = "ancient"

p <- ggplot(psO,aes(x = a2,fill = forcats::fct_rev(psName)))+
  ylab("proportion")+
  geom_bar(stat = "count",position = "fill")+
  scale_fill_manual(values = mypalette2[c(2,4,5,6)],name = "phylostrata")+
  main_theme+
  coord_flip()+
  xlab("abdominal DEG type")+
  scale_x_discrete(expand=c(0,0))+
  guides(fill = guide_legend(reverse=T))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle=-25,hjust=0),
        legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(-0.5,2,0.5,2),"cm"))

psO$a2 = factor(psO$a2,levels = c("conserved queen","conserved worker","non-DE/non-conserved"))
psO$psName[psO$psName=="aculeata"]="hymenoptera"

pH <- ggplot(psO,aes(x = a2,fill = forcats::fct_rev(psName)))+
  ylab("proportion of genes")+
  geom_bar(stat = "count",position = "fill")+
  scale_fill_manual(values = mypalette2[c(2,4,6)],name = "phylostrata")+
  theme(legend.key = element_rect(size=10))+
  main_theme+
  xlab("abdominal caste bias")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  guides(fill = guide_legend(reverse=T))+
  theme(axis.text.x = element_text(angle=-25,hjust=0),
        legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0.5,3.5,0.5,3.5),"cm"))

ggsave(pH,file = "~/GitHub/devnetwork/figures/psCon.png",dpi=300,height=6,width=8)


ggsave(arrangeGrob(pQ,pSoc,p,heights = c(0.4,0.4,0.2)),file="~/GitHub/devnetwork/figures/Fig1.png",height=16,width=10,dpi=350)
ggsave(pSoc,file="~/GitHub/devnetwork/figures/Fig1b.png",height=6,width=10,dpi=300)
ggsave(pQ,file="~/GitHub/devnetwork/figures/Fig1a.png",height=6,width=10,dpi=300)
ggsave(pH,file="~/GitHub/devnetwork/figures/Fig1c.png",height=6,width=6,dpi=300)


#########
###Figure 2a-c
#########
AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

ogg11 = ACUogg
ogg11$abdDE = "non-conserved/non-DE"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "ant worker, bee queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "ant queen, bee worker"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "conserved queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "conserved worker"

FC = merge(AsexRes[[3]][[1]],AcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Mphar",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-conserved/non-DE"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-conserved/non-DE"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-conserved/non-DE"]=0.8

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
  annotate("text",label="A",size=16,x = -8,y=8)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="top",
        legend.text = element_text(size=17),
        legend.title = element_text(size=19,face="bold"))

FC = merge(BsexRes[[3]][[1]],BcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Amel",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-DE/inconsistent"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-DE/inconsistent"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-DE/inconsistent"]=0.8

p2 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  annotate("text",label="B",size=16,x = -8,y=8)+
  ylab("bee caste bias (queen/worker)")+
  xlab("bee sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())

FCb = merge(BsexRes[[3]][[1]],ogg11,by.x="Gene",by.y="gene_Amel")
FC2 = merge(FCb,AsexRes[[3]][[1]],by.x = "gene_Mphar",by.y= "Gene")
FC2$abdDE = factor(FC2$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-conserved/non-DE"))
FC2$alpha = 0.2
FC2$alpha[FC2$abdDE!="non-conserved/non-DE"]=0.8

p3 <- ggplot(FC2,aes(x = -FC.y,y = -FC.x))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  annotate("text",label="C",size=16,x = -8,y=8)+
  ylab("bee sex bias (queen/male)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())


plots <- list(p1,p2,p3)
g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
plots = lapply(plots, function(x) x + theme(legend.position="none"))

p <- arrangeGrob(
  legend,
  do.call(arrangeGrob,c(plots,ncol=3)),
  heights=c(0.1,0.9))

ggsave(p,file="~/GitHub/devnetwork/figures/Fig2a_c.png",width=16,height=6,dpi=300)

p1_blank <- ggplot(FC,aes(x = -FC.x,y = -FC.y,fill=abdDE))+ #Queen-up will be positive
  main_theme+
  scale_fill_manual(values = SexPal,name = "abdominal caste bias")+
  ylab("ant caste bias (queen/worker)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="top",
        legend.text = element_text(size=17),
        legend.title = element_text(size=19,face="bold"))

p2_blank <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  main_theme+
  ylab("bee caste bias (queen/worker)")+
  xlab("bee sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())

p3_blank <- ggplot(FC2,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  ylab("bee sex bias (queen/male)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="none",legend.title = element_blank())

p <- arrangeGrob(p1_blank,p2_blank,p3_blank,nrow=1)

ggsave(p,file="~/GitHub/devnetwork/figures/Fig2a_c_blank.png",width=16,height=6*0.9,dpi=300)

AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

ogg11 = ACUogg
ogg11$abdDE = "non-conserved/non-DE"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "ant worker, bee queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "ant queen, bee worker"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "conserved queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "conserved worker"

FC = merge(AsexRes[[3]][[1]],AcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Mphar",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-conserved/non-DE"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-conserved/non-DE"))
FC = FC[!grepl("ant",FC$abdDE),]
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-conserved/non-DE"]=0.8

p1 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal[c(1,2,5)],name = "abdominal caste bias")+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  scale_alpha_continuous(guide="none")+
  plot2theme+
  ylab("caste bias")+
  xlab("sex bias")+
  ggtitle("ant")+
  ylim(-10,10)+xlim(-10,10)+
  #annotate("text",label="A",size=16,x = -8,y=8)+
  theme(legend.position="top",
        legend.text = element_text(size=17),
        legend.title = element_text(size=19,face="bold"))

FC = merge(BsexRes[[3]][[1]],BcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Amel",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "non-DE/inconsistent"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-DE/inconsistent"))
FC = FC[!grepl("ant",FC$abdDE),]
FC$alpha = 0.2
FC$alpha[FC$abdDE!="non-DE/inconsistent"]=0.8

p2 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal[c(1,2,5)])+
  plot2theme+
  #annotate("text",label="B",size=16,x = -8,y=8)+
  ylab("")+
  xlab("sex bias")+
  ggtitle("bee")+
  ylim(-10,10)+xlim(-10,10)+
  theme(legend.position="none",legend.title = element_blank())


plots <- list(p1,p2)
g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
plots = lapply(plots, function(x) x + theme(legend.position="none"))

p <- arrangeGrob(
  do.call(arrangeGrob,c(plots,ncol=2)),
  legend,
  heights=c(0.9,0.1))

ggsave(p,file="~/GitHub/devnetwork/figures/Fig2a_b.png",width=12,height=6,dpi=300)
########
##Figure 2d
########
DmelSC = merge(sexGenes,ogg11,by="gene_Amel")
DmelSC$abdDE = factor(DmelSC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-DE/inconsistent"))
levels(DmelSC$abdDE)[c(1,2)] = c("conserved queen","conserved worker")
DmelSC$alpha = 0.1
DmelSC$alpha[DmelSC$abdDE!="non-DE/inconsistent"]=0.8

p4 <- ggplot(DmelSC[grepl("conserved",DmelSC$abdDE),],aes(x = abdDE,y=-logFC))+
  #geom_violin(aes(fill=abdDE),alpha=0.8)+
  geom_violin(fill="grey90",trim=FALSE)+
  geom_jitter(width = 0.1,size=0.5,aes(color=abdDE))+
  geom_boxplot(width=0.05,outlier.shape = NA,fill="black",color="black",notch=TRUE,notchwidth = 0.7)+
  plot2theme+
  ylab("fly sex bias")+
  ylim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=0,b=0)))+
  xlab("abdominal caste bias")+
  scale_fill_manual(values=SexPal)+
  scale_color_manual(values=SexPal)+
  theme(legend.position="none")+
  #annotate("text",label="D",size=16,x = 0.75,y=6)+
  stat_summary(geom = "crossbar", width=0.035, fatten=0, size=0.7,color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  coord_cartesian(ylim = c(-8,8))

ggsave(p4,file="~/GitHub/devnetwork/figures/Fig2d.png",width=6,height=6,dpi=300)



plots <- list(p1,p2,p3,p4)
plots[[1]] <- plots[[1]]+theme(legend.margin=margin(t=-900,l=0,r=0,b=0))
g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lwidth <- sum(legend$width)
plots = lapply(plots, function(x) x + theme(legend.position="none",
                                            plot.margin = unit(rep(.75,4),"cm")))

p <- arrangeGrob(
  do.call(arrangeGrob,c(plots,ncol=1)),
  legend,
  widths = unit.c(unit(1, "npc") - lwidth, lwidth))
ggsave(p,file="~/GitHub/devnetwork/figures/Fig2.png",width=8,height=16,dpi=300)


v1 = v2 = rep(0,nrow(DmelSC))
names(v1) = names(v2) = DmelSC$Gene
v1[DmelSC$abdDE=="conserved queen"]=1
v2[DmelSC$abdDE == "conserved worker"]=1
GO_abd <- lapply(list(v1,v2),GOfunc)


GOcb2 <- lapply(GO_abd,function(x) x[c(1:5),c(1,2,6)])
names = c("ant queen, bee worker","ant worker, bee queen","conserved queen","conserved worker")
GOcb3 <- lapply(seq(1,4),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3[c(3,4,1,2)])

t = table(psO$a2,psO$psName)
t = cbind(t,rowSums(t[,-c(1)]))

a = rbind(c(765,93),c(2186,466))
a = rbind(c(259,79),c(2186,466))

aS = merge(Aps,antRes[[2]],by="Gene")
aS = merge(aS,psO,by.x="Gene",by.y="gene_Mphar",all.x=T)
aS$type = "no ortholog"
aS$type[aS$ortholog_found] = "paralog"
aS$type[aS$OGG_found] = "ortholog, non-conserved"
aS$type[grepl("conserved",aS$abdDE)] = "ortholog, conserved"
ggplot(aS,aes(x = type,fill=psName.x))+
  geom_bar(stat="count",position="fill")+
  facet_grid(. ~ abdomen)

a = rbind(c(1020,1206-1020),c(7825-1020,10011-7825-1206+1020))


psO = merge(AllPS_sum,ogg11,by="OGGacu")

aS = merge(Bps,beeRes[[2]],by="Gene")
aS = merge(aS,psO,by.x="Gene",by.y="gene_Amel",all.x=T)
aS$type = "no ortholog"
aS$type[aS$ortholog_found] = "paralog"
aS$type[aS$OGG_found] = "ortholog, non-conserved"
aS$type[grepl("conserved",aS$abdDE)] = "ortholog, conserved"
ggplot(aS,aes(x = type,fill=psName.x))+
  geom_bar(stat="count",position="fill")


Aps$species = "M. pharaonis"
Bps$species = "A. mellifera"
Aps_FC = merge(Aps,antRes[[1]],by="Gene")
Bps_FC = merge(Bps,beeRes[[1]],by="Gene")

ggplot(Bps_FC,aes(x = psName,fill=larva))+
  geom_bar(stat = "count",position="fill")

for (i in 7:11){
  Aps_FC[,i] = Aps_FC[,i] - median(Aps_FC[,i])
  Bps_FC[,i] = Bps_FC[,i] - median(Bps_FC[,i])
}
Apm = melt(Bps_FC,id.vars = c("Gene","ODBgene","OGGacu","ps","psName","species"))

p5 <- ggplot(Aps_FC,aes(x = psName,y=-abdomen,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape = NA)+
  coord_cartesian(ylim = c(-6,6))+
  main_theme+
  scale_fill_manual(values = rev(mypalette),name = "phylostrata")+
  ylab("log2 fold-change\n(queen abdomen/worker abdomen)")+
  geom_hline(yintercept =0,linetype="dashed")+
  theme(legend.position=c(0.8,0.8),
        #legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank())

p6 <- ggplot(Bps_FC[!is.na(Bps_FC$psName),],aes(x = psName,y=-abdomen,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape = NA)+
  coord_cartesian(ylim = c(-6,6))+
  main_theme+
  scale_fill_manual(values = rev(mypalette),name = "phylostrata")+
  ylab("log2 fold-change\n(queen abdomen/worker abdomen)")+
  geom_hline(yintercept =0,linetype="dashed")+
  theme(legend.position=c(0.8,0.8),
        #legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank())

png("~/GitHub/devnetwork/figures/Fig2.png",width=3000,height=6000,res=300)
grid.arrange(grid_arrange_shared_legend(p1,p2,p3,p4,position = "right",ncol=1))
dev.off()

Bps_FC = Bps_FC[!is.na(Bps_FC$psName),]
binom.test(sum(Bps_FC$abdomen[Bps_FC$psName=="novel"] < 0),sum(Bps_FC$psName=="novel"),0.5)

AllP = rbind(Aps_FC,Bps_FC)
AllP$species = factor(AllP$species,levels=c("M. pharaonis","A. mellifera"))
AllP=droplevels(AllP[!is.na(AllP$psName),])
AllP$psName=as.character(AllP$psName)
AllP$psName[AllP$psName=="bee"|AllP$psName=="ant"]="ant/bee"
AllP$psName=factor(AllP$psName,levels=c("old","insect","hymenoptera","aculeata","ant/bee","novel"))

a = rev(mypalette)

p5 <- ggplot(AllP,aes(x = species,y=-abdomen,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape = NA)+
  coord_cartesian(ylim = c(-6,6))+
  main_theme+
  scale_fill_manual(values = a,name = "phylostrata")+
  ylab("log2 fold-change\n(queen abdomen/worker abdomen)")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position='right',
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        panel.border=element_blank(),
        axis.line.x=element_line("black",size=1),
        axis.line.y=element_line("black",size=1),
        legend.margin = margin(t=290,r=0,l=-50,b=0))


g <- ggplotGrob(p1 + theme(legend.position = "top"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

pAll = lapply(list(p1,p2,p3,p4), function(x) x+ theme(legend.position="none"))
pAll = c(pAll,legend)

png("~/GitHub/devnetwork/figures/Fig2all.png",width=3000,height=2000,res=300)
grid.arrange(
  pAll[[1]],legend,
  nrow=2
)
dev.off()

#########
##Figure 3
#########
antQG <- idQG(antPl,7,factorA,antT)
beeQG <- idQG(beePl,7,factorB,beeT)
aG <- commonGenes(antQG,antT)
bG <- commonGenes(beeQG,beeT)
antG = aG$Gene[aG$KeepNum > 0.6]
beeG = bG$Gene[bG$KeepNum > 0.6]
ogg11$antBic = ogg11$beeBic = 0
ogg11$antBic[ogg11$gene_Mphar %in% antG]=1
ogg11$beeBic[ogg11$gene_Amel %in% beeG] = 1
ogg11$totBic = ogg11$antBic + ogg11$beeBic

#Make heatmaps of an example bicluster
clusters <-as.data.frame(beePl[[1]]@NumberxCol)
colnames(clusters) = colnames(beeT)
genes = as.data.frame(beePl[[1]]@RowxNumber)
rownames(genes) = rownames(beeT)

sizes = colSums(beePl[[1]]@RowxNumber)
clusters = clusters[sizes >= 25,]
genes = genes[,sizes >= 25]

hc = hclust(dist(t(clusters)))
dend <- as.dendrogram(hc)
#order so QG is on the right
dend[[2]][[2]] = list(dend[[2]][[2]][[2]],dend[[2]][[2]][[1]])
dend[[2]][[2]][[2]] = list(dend[[2]][[2]][[2]][[2]],dend[[2]][[2]][[2]][[1]])
hc$order = order.dendrogram(dend)
clusters = clusters[,order.dendrogram(dend)]

#Normalize quantiles for expression data
#From https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


gK = genes[rowSums(genes) == 1,]
geneCluster = apply(gK,1,function(x) seq(1,length(x))[x])
gK = rownames(gK)[order(geneCluster,decreasing=TRUE)]

m1 = beeT[gK,order.dendrogram(dend)]
m1 = log(m1+sqrt(m1^2+1))
m1 = quantile_normalisation(m1)
m1 = m1 - rowSums(m1)/ncol(m1)
m1 = m1 - colSums(m1)/nrow(m1)
d = melt(as.matrix(m1))

## set color representation for specific values of the data distribution
quantile_range <- quantile(d$value, probs = seq(0,1,by=0.1))
color_palette <- colorRampPalette(c("slateblue4","white","firebrick1"))(length(quantile_range) - 1)
d$value_disc <- findInterval(d$value, quantile_range, all.inside = TRUE)
label_text <- rollapply(seq(0,100,by=10), width = 2, by = 1, FUN = function(i) paste(i, collapse = " - "))

pBic <- ggplot(d,aes(x = Var1,y=Var2,fill=factor(value_disc)))+
  geom_raster()+
  main_theme+
  scale_fill_manual(values=color_palette,name = "expression\npercentile",
                    labels = label_text,guide=guide_legend(title.position = "top",title.hjust = 0.5))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=7,face="bold"),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2,"cm"),
        legend.margin = margin(t=150,l=0,r=0,b=0))

png("~/Github/devnetwork/figures/Fig3a.png",height=2000,width=2000,res=600)
p1
dev.off()

#########
##Figure 3b
#########
aG2 = merge(aG,ogg11,by.x="Gene",by.y="gene_Mphar",all.x=TRUE)
aG2$OGGacu = as.character(aG2$OGGacu)
aG2$OGGacu[is.na(aG2$OGGacu)] = as.character(aG2$Gene)[is.na(aG2$OGGacu)]
antG2 = aG2$OGGacu[aG2$KeepNum > 0.5]
bG2 = merge(bG,ogg11,by.x="Gene",by.y="gene_Amel",all.x=TRUE)
bG2$OGGacu = as.character(bG2$OGGacu)
bG2$OGGacu[is.na(bG2$OGGacu)] = as.character(bG2$Gene)[is.na(bG2$OGGacu)]
beeG2 = bG2$OGGacu[bG2$KeepNum > 0.5]

p2 <- venn.diagram(list("M. pharaonis"=antG2,"A. mellifera"=beeG2,"1-1 orthologs"=ogg11$OGGacu),filename=NULL,
             fill = c("red","blue","white"),
             imagetype = "png",
             cex=1.5,
             cat.cex=1.6,
             height=3000,
             width=6000,
             cat.pos = c(-20,20,0),
             fontface="bold",
             cast.dist = c(100,0,30),
             cat.fontface = c("italic","italic","plain"),
             margin = 0.05)


p2b <- venn.diagram(list("M. pharaonis"=antG2[antG2 %in% ogg11$OGGacu],"A. mellifera"=beeG2[beeG2 %in% ogg11$OGGacu]),filename=NULL,
                   imagetype = "png",
                   main = "number of 1-1 orthologs in bicluster",
                   cex=1.5,
                   main.cex = 2.2,
                   main.pos = c(0.5,0.15),
                   cat.cex=1.6,
                   height=3000,
                   width=7000,
                   cat.pos = c(-15,15),
                   fontface="bold",
                   cat.dist = rep(0.04,2),
                   cat.fontface = c("italic","italic"),
                   margin = 0.05)

cats = table(ogg11$abdDE[ogg11$totBic==2])
slices <- cats[c(3,4,2,1,5)]
n = names(cats)[c(3,4,2,1,5)]
dn = data.frame(slices)
dn$x="yes"

p3 <- ggplot(dn,aes(fill=Var1,y=Freq,x=x))+
  geom_bar(stat="identity")+
  main_theme+
  scale_fill_manual(values=SexPal,name="caste bias")+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_text(face="bold"),
        legend.margin = margin(t=200,l=0,r=0,b=0))

png("~/GitHub/devnetwork/figures/Fig3b2.png",height=2000,width=1500,res=300)
p3
dev.off()

########
##Figure 3c,d
#######
getConn <- function(tpm,genes,DEtests,DEres){
  tpm = tpm[rownames(tpm) %in% genes,]
  adj1 = cor(t(tpm))^6
  conns = data.frame(Gene=rownames(tpm),conn=rowSums(adj1,na.rm = T))
  bDG =DEres[[1]]
  bDG$DEcat = DEres[[2]]$abdomen
  bFDR = data.frame(Gene=rownames(DEtests[["8_gaster"]]),FDR = DEtests[["8_gaster"]]$FDR)
  bDG = merge(bDG,bFDR,by="Gene")
  beeGr = bDG[bDG$Gene %in% genes,]
  beeGrC = merge(beeGr,conns,by="Gene")
  beeGrC$DEcat = factor(beeGrC$DEcat,levels = c("queen","worker","nonDE"))
  levels(beeGrC$DEcat) = c("queen","worker","non-DE")
  return(beeGrC)
}

beeGrC = getConn(beeT,beeG,beeTests,beeRes)
antGrC = getConn(antT,antG,antTests,antRes)


antGrC = merge(antGrC,Aps,by="Gene")
beeGrC = merge(beeGrC,Bps,by="Gene")

beeGrC = merge(beeGrC,aName,by="Gene")
beeGrC$size=2
beeGrC$size[grepl("Smaug|vitellogenin|vasa|ovo",beeGrC$GeneName)]=4
b = beeGrC[beeGrC$conn > quantile(beeGrC$conn,0.9) & beeGrC$DEcat == "queen" & beeGrC$abdomen < -2 & beeGrC$Gene %in% ogg11$gene_Amel[ogg11$abdDE=="conserved queen"],]

ann <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/MpharAnn.csv")
antGrC = merge(antGrC,ann[,c(1,2)],by="Gene")
colnames(antGrC)[ncol(antGrC)]="GeneName"
antGrC$size=2
antGrC$size[grepl("Smaug|Vitellogenin-2|vasa|Vitellogenin receptor",antGrC$GeneName)]=4
b = beeGrC[beeGrC$conn > quantile(beeGrC$conn,0.9) & beeGrC$DEcat == "queen" & beeGrC$abdomen < -2 & beeGrC$Gene %in% ogg11$gene_Amel[ogg11$abdDE=="conserved queen"],]
genes = b$SwissProt
antGrC$conn = antGrC$conn/max(antGrC$conn)
beeGrC$conn = beeGrC$conn/max(beeGrC$conn)

allC = rbind(antGrC,beeGrC)
levels(allC$psName)[4:7] = "aculeata &\nyounger"
levels(allC$psName)[1] = "ancient"
allC = merge(allC,ogg11,by = "OGGacu",all.x=T)
allC$abdDE = factor(allC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","non-conserved/non-DE"))
allC$abdDE[is.na(allC$abdDE)] = "non-conserved/non-DE"
allC$alpha = 0.9
allC$alpha[allC$abdDE=="non-conserved/non-DE"]=0.8
allC$species=factor(allC$species,levels = c("M. pharaonis","A. mellifera"))
levels(allC$species) = c("ant","bee")

p1 <- ggplot(allC,aes(x = conn,y=-abdomen))+
  geom_point(aes(fill = DEcat),pch=21,color="black")+
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
        strip.background = element_rect(color=NA,fill=NA))

ggsave(p1,file = "~/GitHub/devnetwork/figures/bicConn1.png",height=6,width=10,dpi=300)


p1 <- ggplot(allC,aes(x = conn,y=-abdomen))+
  geom_point(aes(fill = DEcat,size=size),pch=21,color="black")+
  scale_fill_manual(values = SexPal[c(1,2,5)],name="DEG type")+
  scale_alpha_continuous(guide="none")+
  scale_size_continuous(guide="none")+
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
        strip.background = element_rect(color=NA,fill=NA))

ggsave(p1,file = "~/GitHub/devnetwork/figures/bicConn.png",height=6,width=10,dpi=300)

p2 <- ggplot(allC,aes(x = psName,y=conn))+
  geom_violin(aes(fill=psName),alpha = 0.9,trim=FALSE)+
  ylab("connectivity within bicluster")+
  geom_boxplot(width=0.15,outlier.shape = NA,fill="black")+
  scale_fill_manual(values=rev(mypalette2[c(2,4,5,6)]))+
  xlab("phylostrata")+
  facet_grid(. ~ species)+
  main_theme+
  scale_y_log10(breaks=c(1,10,100))+
  stat_summary(geom = "crossbar", width=0.05, fatten=0, size=1,color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  stat_n_text(y.pos = -0.9,size = 5)+
  theme(legend.position = "none",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        plot.margin = unit(rep(1,4),"cm"),
        axis.text.x = element_text(angle=-25,hjust=0),
        strip.background = element_rect(color="black",fill="darkgrey"))

p <- arrangeGrob(p1,p2,ncol=1)
pBic <- pBic+theme(legend.position = "right",
                 legend.margin = margin(t=5,r=5,l=5,b=5),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=15),
                 legend.key.size = unit(0.65,"cm"),
                 plot.margin = margin(t=15,l=15,r=15,b=15),
                 axis.title = element_text(size=20,face="bold"))+
  guides(fill = guide_legend(reverse=T))+
  ylab("samples")+
  xlab("genes")

pL = arrangeGrob(pBic,p,nrow=2,heights = c(0.3,0.7))

ggsave(pL,file="~/GitHub/devnetwork/figures/Fig3.png",width=12,height=16,dpi=300)


p1_blank <- ggplot(allC,aes(x = conn,y=-abdomen))+
  ylab("abdominal caste bias (queen/worker)")+
  xlab("connectivity in abdomen bicluster")+
  facet_grid(. ~ species)+
  main_theme+
  scale_x_log10()+
  theme(legend.position = c(0.15,0.8),
        legend.text = element_text(size=15),
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        plot.margin = unit(rep(1,4),"cm"),
        strip.background = element_rect(color="black",fill="darkgrey"))


p2_blank <- ggplot(allC,aes(x = psName,y=conn))+
  ylab("connectivity within bicluster")+
  xlab("phylostrata")+
  facet_grid(. ~ species)+
  main_theme+
  scale_y_log10(breaks=c(1,10,100))+
  theme(legend.position = "none",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        plot.margin = unit(rep(1,4),"cm"),
        axis.text.x = element_text(angle=-25,hjust=0),
        strip.background = element_rect(color="black",fill="darkgrey"))

ggsave(p1_blank,file="~/GitHub/devnetwork/figures/Fig3b_blank.png",width=12,height=16*0.35,dpi=300)
ggsave(p2_blank,file="~/GitHub/devnetwork/figures/Fig3c_blank.png",width=12,height=16*0.35,dpi=300)
ggsave(p1,file="~/GitHub/devnetwork/figures/Fig3b.png",width=12,height=16*0.35,dpi=300)
ggsave(p2,file="~/GitHub/devnetwork/figures/Fig3c.png",width=12,height=16*0.35,dpi=300)


wilcox.test(allC$conn[allC$psName=="ancient" & allC$species == "A. mellifera"],allC$conn[allC$psName=="aculeata &\nyounger"& allC$species == "A. mellifera"])


png("~/GitHub/devnetwork/figures/Fig3all.png",width=4000,height=4000,res=300)
grid.arrange(p1+theme(legend.position = c(0.15,0.35),
                      legend.margin = margin(t=5,r=5,l=5,b=5),
                      legend.text = element_text(size=12),
                      legend.title = element_text(size=15),
                      legend.key.size = unit(0.65,"cm"),
                      plot.margin = margin(t=15,l=15,r=15,b=15),
                      axis.title = element_text(size=20,face="bold"))+
               guides(fill = guide_legend(reverse=T))+
               ylab("samples")+
               xlab("genes"),
             p4+
               theme(legend.position = c(0.8,0.2),
                     axis.title = element_text(size=20,face="bold"),
                     legend.background = element_blank())+
               theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=0,b=0)),
                     legend.title = element_text(face = "bold",size=17))+
               theme(plot.margin=margin(t=15,r=15,l=15,b=15)),
             grobTree(p2b),
             p5+
               theme(legend.position = c(0.8,0.2),
                     legend.background = element_blank(),
                     axis.title = element_text(size=20,face="bold"))+
               theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=0,b=0)),
                     legend.title = element_text(face = "bold",size=17))+
               theme(plot.margin=margin(t=15,r=15,l=15,b=15)),
             ncol=2)
grid.rect(x=0.172,y=0.968,width=0.076,height=0.032,gp=gpar(color="black",lwd=5,fill=NA))
grid.text(label="A",x=0.07,y=0.95,gp=gpar(cex=3))
grid.text(label="C",x=0.62,y=0.95,gp=gpar(cex=3))
grid.text(label="B",x=0.07,y=0.48,gp=gpar(cex=3))
grid.text(label="D",x=0.62,y=0.45,gp=gpar(cex=3))
dev.off()


a = rbind(c(463+1159,705),c(5145+873-1159-463,nrow(ant) - (705+1159+463+5145+873)))
a = rbind(c(463+873,449),c(5145-873+1159-463,nrow(bee) - (705+1159+463+5145+873)))

a = rbind(c(463,873),c(1159,5145))

BpsC = merge(beeGrC,Bps,by="Gene")
ApsC = merge(antGrC,Aps,by="Gene")
BpsC$psName[grepl("aculeata|bee|novel",BpsC$psName)] = "bee"

ggplot(BpsC,aes(x = psName,y=conn))+
  geom_boxplot(notch=TRUE)

ApsC$psName[grepl("aculeata|ant|novel",ApsC$psName)] = "ant"

ggplot(ApsC,aes(x = psName,y=conn))+
  geom_boxplot(notch=TRUE)

#Find some interesting genes
res = merge(antRes[[1]],ogg11,by.x = "Gene",by.y = "gene_Mphar")
res = merge(res,beeRes[[1]],by.x = "gene_Amel",by.y = "Gene")
res = res[,c(1,2,7,9:11,18)]
colnames(res) = c("gene_Amel","gene_Mphar","antFCabd","caste bias","beeAbdBic","antAbdBic","beeFCabd")
res = res[,c(2,1,3,7,5,6,4)]
res$antFCabd=-res$antFCabd
res$beeFCabd=-res$beeFCabd
res = merge(res,ann,by.x="gene_Mphar",by.y = "gene")[,c(1:7,17,30)]
res = res[res$DescriptionSP!="-",]
res=res[order(res$beeFCabd^2+res$antFCabd^2,decreasing=TRUE),]
write.csv(res,file="~/GitHub/devnetwork/figures/conservedcaste.csv")

#########
##Figure 4a
#########
antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]])
antSB = euclDist(antSocRes[[1]])
beeSB = euclDist(beeSocRes[[1]])
antCB_pval = geomMean(antRes_allstage[[3]])
beeCB_pval = geomMean(beeRes_allstage[[3]])
antSB_pval = geomMean(antSocRes[[3]])
beeSB_pval = geomMean(beeSocRes[[3]])

antBias = merge(antCB,antSB,by="Gene")
beeBias = merge(beeCB,beeSB,by="Gene")
antCB$type=beeCB$type="caste"
antSB$type=beeSB$type="behavior"
antSB$cb_noAdult=antSB$cb
beeSB$cb_noAdult=beeSB$cb
antBias2 = rbind(antCB,antSB)
beeBias2 = rbind(beeCB,beeSB)
cbBias = merge(beeCB,ACUogg,by.x="Gene",by.y="gene_Amel")
cbBias = merge(cbBias,antCB,by.x="gene_Mphar",by.y="Gene")
cbBias = cbBias[!is.na(cbBias$OGGacu),]
socBias = merge(beeSB,ACUogg,by.x="Gene",by.y="gene_Amel")
socBias = merge(socBias,antSB,by.x="gene_Mphar",by.y="Gene")
socBias = socBias[!is.na(socBias$OGGacu),]

antCor = cor(t(antT[rowSums(antT) > 0,]))^6
beeCor = cor(t(beeT[rowSums(beeT) > 0,]))^6
antConn = data.frame(Gene = rownames(antCor),kTotal = rowSums(antCor))
beeConn = data.frame(Gene = rownames(beeCor), kTotal = rowSums(beeCor))

##########
##Fig 4e
##########
cbAps = merge(antBias2,Aps,by="Gene")
cbBps = merge(beeBias2,Bps,by="Gene")
cbAps = merge(antConn,cbAps,by="Gene")
cbBps = merge(beeConn,cbBps,by="Gene")

cbBps$kTotal = cbBps$kTotal/max(cbBps$kTotal)
cbAps$kTotal = cbAps$kTotal/max(cbAps$kTotal)

levels(cbAps$psName)[1] = levels(cbBps$psName)[1]= "ancient"
pA <- ggplot(cbAps[cbAps$type=="caste",],aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("ant")+
  xlab("")+
  coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")


pB <- ggplot(cbBps[cbBps$type=="caste" & !is.na(cbBps$psName),],aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  plot2theme+
  ggtitle("bee")+
  xlab("")+
  coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("")+
  theme(axis.text.x = element_text(angle = -25,hjust=0),
        legend.position="none")

p <- arrangeGrob(arrangeGrob(pA,pB,ncol=2),
                 textGrob("estimated evolutionary age",x = 0.5,y=1.5,gp = gpar(cex=1.8,fontface="bold")),
                 ncol=1,heights = c(0.95,0.05))

ggsave(p,file = "~/GitHub/devnetwork/figures/cb_phylo.png",height=6,width=10,dpi=300)

cbP = rbind(cbAps,cbBps)

cbP$species = factor(cbP$species,levels=c("M. pharaonis","A. mellifera"))
cbP=droplevels(cbP[!is.na(cbP$psName),])
cbP$psName=as.character(cbP$psName)
cbP$psName[cbP$psName=="bee"|cbP$psName=="ant"]="ant/bee"
cbP$psName=factor(cbP$psName,levels=c("old","insect","hymenoptera","aculeata","ant/bee","novel"))
cbP$type = factor(cbP$type,levels = c("caste","behavior"))
levels(cbP$psName)[1] = "ancient"

p1ph <- ggplot(cbP[cbP$type=="caste",],aes(x = species,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position="right",
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


p1F <- ggplot(cbP[cbP$type=="caste",],aes(x = psName,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  facet_grid(. ~ species)+
  coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  xlab("phylostata")+
  theme(strip.text = element_text(size=20,face="bold.italic"),
    axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
    legend.position = "none",
    plot.margin = unit(rep(1,4),"cm"),
    axis.text.x = element_text(angle=-25,hjust=0),
    strip.background = element_rect(color="black",fill="darkgrey"))

ggsave(p1F,filename = "~/GitHub/devnetwork/figures/Fig4ef.png",height=8,width=12,dpi=300)


p1ph <- ggplot(cbP[cbP$type=="caste",],aes(x = species,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,2))+
  scale_y_continuous(breaks = c(0,1,2,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position="right",
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


p2ph <- ggplot(cbP[cbP$type=="behavior",],aes(x = species,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,3))+
  scale_fill_manual(values = rev(mypalette2),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall behavior bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position="right",
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        legend.margin = margin(t=460,r=0,l=-20,b=0),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


plots <- list(p1ph,p2ph)
g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
plots = lapply(plots, function(x) x + theme(legend.position="none"))

p <- arrangeGrob(
  legend,
  do.call(arrangeGrob,c(plots,ncol=2)),
  heights=c(0.15,0.85),
  ncol=1)

ggsave(p,file = "~/GitHub/devnetwork/figures/Fig4ef.png",height=6,width=10,dpi=300)
cbAps$kTotal = cbAps$kTotal/max(cbAps$kTotal)
cbBps$kTotal = cbBps$kTotal/max(cbBps$kTotal)

p1C <- ggplot(cbAps[cbAps$type=="caste",],aes(x = kTotal,y=cb_noAbd))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("scaled network connectivity")+
  ggtitle("ant")+
  ylab("overall caste bias")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  theme(legend.position="none")

p2C <- ggplot(cbAps[cbAps$type=="behavior",],aes(x = kTotal,y=cb))+
  geom_hex()+
  main_theme+
  xlab("scaled network connectivity")+
  ylab("overall ant behavior bias")+
  scale_x_log10()+
  ggtitle("bee")+
  geom_smooth(se=FALSE,color="black")+
  theme(legend.position="none",
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        legend.margin = margin(t=460,r=0,l=-20,b=0),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


p3C <- ggplot(cbBps[cbBps$type=="caste",],aes(x = kTotal,y=cb_noAbd))+
  geom_hex(bins=70)+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  scale_fill_gradient(low = "blue",high="red")+
  xlab("scaled network connectivity")+
  ylab("")+
  ggtitle("bee")+
  scale_x_log10(breaks = c(0.01,0.1,1))+
  theme(legend.position="none")

p4C <- ggplot(cbBps[cbBps$type=="behavior",],aes(x = kTotal,y=cb))+
  geom_hex()+
  main_theme+
  xlab("total network connectivity")+
  ylab("overall bee behavior bias")+
  scale_x_log10()+
  theme(legend.position="none",
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
  
p <- arrangeGrob(p1C,p3C,ncol=2)
ggsave(p,file="~/GitHub/devnetwork/figures/Fig4ac.png",height=6,width=12,dpi=300)

cor.test(cbBps$kTotal[cbBps$type=="caste"],cbBps$cb_noAbd[cbBps$type=="caste"],method="spearman")
cor.test(cbAps$kTotal[cbAps$type=="caste"],cbAps$cb_noAbd[cbAps$type=="caste"],method="spearman")



p <- arrangeGrob(arrangeGrob(p1C,p2C,p3C,p4C,ncol=1),arrangeGrob(p1ph,p2ph,ncol=1),ncol=2)
ggsave(p,file = "~/GitHub/devnetwork/figures/Fig4.png",height=12,width=10,dpi=300)

png("~/GitHub/devnetwork/figures/Fig4e2.png",width=3500,height=2000,res=300)
grid.arrange(grid_arrange_shared_legend(p1,p2,position="right"))
dev.off()

png("~/GitHub/devnetwork/figures/Fig4e_horiz.png",width=3500,height=1750,res=300)
grid.arrange(p1,p2,ncol=2)
dev.off()

p1 <- ggplot(antBias,aes(x=cb_noAdult.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  main_theme+xlab("ant caste bias")+ylab("ant behavior bias")

p2 <- ggplot(beeBias,aes(x=cb_noAdult.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  main_theme+xlab("bee caste bias")+ylab("bee behavior bias")

p3 <- ggplot(cbBias,aes(x=cb_noAdult.x,y=cb_noAdult.y))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  main_theme+xlab("bee caste bias")+ylab("ant caste bias")

p4 <- ggplot(beeBias,aes(x=cb.x,y=cb.y))+
  geom_point(size=2,alpha=0.3)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  scale_x_continuous(breaks = c(0,1,2))+
  main_theme+xlab("bee behavioral bias")+ylab("ant behavior bias")

png("~/GitHub/devnetwork/figures/Fig4a_d.png",height=3000,width=3000,res=300)
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()




##########
###Fig 5a-d (or supplement)
##########
bT = merge(tau,beeBias2,by="Gene")
bT$species = "bee"
bT$type = factor(bT$type,levels = c("caste","behavior"))
levels(bT$type) = c("caste bias","behavior bias")

p1 <- ggplot(bT,aes(x=tau,y=cb))+
  geom_point(size=1,alpha=0.6,pch=21,aes(fill=type),color="grey60")+
  geom_smooth(se=FALSE,aes(color=type))+
  scale_color_manual(values = SexPal)+
  scale_fill_manual(values = SexPal)+
  main_theme+xlab("tissue specificity")+ylab("bee bias")+
  theme(legend.position=c(0.2,0.6),
        legend.text = element_text(size=20),
        legend.title=element_blank())

bT2 = merge(bT,ACUogg,by.x="Gene",by.y="gene_Amel")
bT3 = merge(bT2,antBias2,by.x=c("gene_Mphar","type"),by.y=c("Gene","type"))
bT3$species = "ant"

bT3$type = factor(bT3$type,levels = c("caste","behavior"))
levels(bT3$type) = c("caste bias","behavior bias")

p2 <- ggplot(bT3,aes(x=tau,y=cb.y))+
  geom_point(size=1,alpha=0.6,pch=21,aes(fill=type),color="grey60")+
  geom_smooth(se=FALSE,aes(color=type))+
  scale_color_manual(values = SexPal)+
  scale_fill_manual(values = SexPal)+
  main_theme+xlab("tissue specificity")+ylab("ant bias")+
  theme(legend.position="none")

png("~/GitHub/devnetwork/figures/Fig5a_b.png",width=4000,height=2000,res=300)
grid.arrange(p1,p2,ncol=2)
dev.off()

F4b <- ggplot(bT,aes(x=tau,y=cb))+
  geom_point(size=1,alpha=0.6,pch=21,aes(fill=type),color="grey60")+
  geom_smooth(se=FALSE,aes(color=type))+
  scale_color_manual(values = specPal,name="bias type")+
  scale_fill_manual(name = "bias type",values = specPal)+
  main_theme+
  xlab(expression(bold("tissue specificity"~(tau))))+
  ylab("overall bee bias")+
  theme(legend.position=c(0.2,0.7),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        plot.margin = unit(c(.75,.75,.75,.75),"cm"))

F4a <- ggplot(bT3,aes(x=tau,y=cb.y))+
  geom_point(size=1,alpha=0.6,pch=21,aes(fill=type),color="grey60")+
  geom_smooth(se=FALSE,aes(color=type))+
  scale_color_manual(values = specPal,name="bias type")+
  scale_fill_manual(values = specPal,name="bias type")+
  main_theme+xlab(expression(bold("tissue specificity"~(tau))))+ylab("overall ant bias")+
  theme(legend.position=c(0.2,0.7),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        plot.margin = unit(c(.75,.75,.75,.75),"cm"))

png("~/GitHub/devnetwork/figures/Fig4a_b.png",width=4000,height=2000,res=300)
grid.arrange(F4a,F4b,ncol=2)
dev.off()

p1ph <- ggplot(cbAps,aes(x = type,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,3))+
  xlab("bias type")+
  scale_fill_manual(values = rev(mypalette),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall ant bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position=c(0.8,0.8),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        plot.margin = unit(c(.75,.75,.75,.75),"cm"))

p2ph <- ggplot(cbBps[!is.na(cbBps$psName),],aes(x = type,y=cb,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,3))+
  xlab("bias type")+
  scale_fill_manual(values = rev(mypalette),name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall bee bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position=c(0.8,0.8),
        #legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        #legend.background = element_blank(),
        #legend.margin = margin(t=220,r=0,l=0,b=0), #460 for vertical
        plot.margin = unit(c(.75,.75,.75,.75),"cm"))

png("~/GitHub/devnetwork/figures/Fig4.png",width=4000,height=4000,res=325)
grid.arrange(F4a,p1ph,F4b,p2ph,ncol=2)
grid.text("A",x=0.1,y=0.93,gp=gpar(cex=3,fontface="bold"))
grid.text("B",x=0.1,y=0.43,gp=gpar(cex=3,fontface="bold"))
grid.text("C",x=0.6,y=0.93,gp=gpar(cex=3,fontface="bold"))
grid.text("D",x=0.6,y=0.43,gp=gpar(cex=3,fontface="bold"))
dev.off()




antBias2$devel = beeBias2$devel = 0
antBias2$devel[antBias2$Gene %in% antDevel] = 1
beeBias2$devel[beeBias2$Gene %in% beeDevel] = 1


#Calculate a gene-wise "developmental index", the euclidean distance of log fold-change between each developmental stage
develIndex <- function(counts,factor){
  f = droplevels(factor[factor$tissue=="larva"|factor$tissue=="egg",])
  i = 1
  results = list()
  f$stage = as.integer(as.character(f$stage))
  for (lev1 in c(1:6)){
    for (lev2 in c(1:6)){
      if (lev2 <= lev1){
        next;
      }
      f1 = droplevels(f[f$stage==lev1|f$stage==lev2,])
      c = counts[,colnames(counts) %in% f1$sample]
      design <- model.matrix(~stage+colony,data = f1)
      res <- EdgeR(c,design,2) 
      res$Gene = rownames(res)
      res = res[,c("logFC","Gene")]
      results[[i]] = res
      i = i+1
    }
  }
  allRes <- join_all(results,by="Gene",type='inner')
  devIndex <- apply(allRes[,-c(2)],1,function(x) sqrt(sum(x^2)))
  names(devIndex) = allRes$Gene
  return(devIndex)
}

devFCBee <- develIndex(bee,factorB)
devFCAnt <- develIndex(ant,factorA)
B1 = data.frame

ggplot(beeBias2,aes(x=type,y=cb,fill=as.factor(devel)))+
  geom_boxplot(notch=TRUE)+
  main_theme+
  xlab("comparison")+
  ylab("overall bias")





tissueExpr <- beeTissue[[2]]
tissueExpr$Gene = rownames(tissueExpr)
topTissue = apply(tissueExpr[,c(1:12)],1,function(x) colnames(tissueExpr)[1:12][x==max(x)])
topTissue <- lapply(topTissue,function(x){
  if (length(x) > 1) NA
  else x
})
t = unlist(topTissue)
tissueExpr$topTissue = t
dBTissue <- merge(bT,tissueExpr[,c(13,14)],by="Gene")
dBTissue = dBTissue[!is.na(dBTissue$tau),]
dBTissue_topTau <- dBTissue[dBTissue$tau > quantile(dBTissue$tau,0.9),]
dBTissue_topTau$gland_type = "conserved"
dBTissue_topTau$gland_type[dBTissue_topTau$topTissue=="hypophar" | dBTissue_topTau$topTissue== "mandibular" |
                             dBTissue_topTau$topTissue == "nasonov" | dBTissue_topTau$topTissue == "sting" |
                             dBTissue_topTau$topTissue == "antenna"] = "novel"
dBTissue_topTau$gland_type = as.factor(dBTissue_topTau$gland_type)

png("~/GitHub/devnetwork/figures/Fig5c.png",height = 2000, width = 2000, res = 300)
ggplot(dBTissue_topTau,aes(x = type,y=cb_noAdult,fill = gland_type))+
  geom_boxplot(notch=TRUE,outlier.shape = NA)+
  main_theme+
  scale_fill_manual(values = rev(mypalette[c(1,6)]),name = "gland type")+
  ylab("bee bias")+
  xlab("")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = c(0.8,0.8),
        legend.title = element_text(size=17,face="bold"),
        legend.text = element_text(size = 15))+
  coord_cartesian(ylim = c(0,3))
dev.off()

dBTissue_topTau$topTissue = factor(dBTissue_topTau$topTissue,levels = c("antenna","hypophar","mandibular","nasonov","sting",
                                                                        "abdomen","brain","digestive","malpighian","midgut","muscle","thor_gang"))

levels(dBTissue_topTau$topTissue)[12] = "ganglion"

png("~/GitHub/devnetwork/figures/topTissue_cb.png",height = 2000, width = 3000, res = 300)
ggplot(dBTissue_topTau,aes(x = topTissue,y=cb_noAdult,fill = gland_type))+
  geom_boxplot(outlier.shape=NA)+
  facet_grid(. ~ type)+
  main_theme+
  ylab("bee bias")+
  coord_cartesian(ylim = c(0,4))+
  scale_fill_manual(values = rev(mypalette[c(1,6)]),name="gland type")+
  xlab("tissue with highest expression")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = c(0.85,0.8),
        legend.title = element_text(size=17,face="bold"),
        legend.text = element_text(size = 15),
        strip.text = element_text(size=15,face="italic"),
        plot.margin = unit(c(0.5,1,0.5,0.5),"cm"))
dev.off()

######
##GO analysis
#####
go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
library(topGO)
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

selectDE <- function(score){
  return(score == 1)
}


selectConn <- function(score){
  return(score > quantile(score,0.9))
}

#GSEA analysis. Takes in vector of all genes, with a score for social connection strength
GSEAfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectConn,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="decreasing")
  allRes <- GenTable(GOdata,KS=resultKS)
  return(allRes) 
}

GOfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectDE,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata,KS=resultKS)
  return(allRes) 
}

v1 = antCB$cb
v2 = antSB$cb
names(v1) = names(v2) = antCB$Gene
# v3 = beeCB$cb
# v4 = beeSB$cb
# names(v3) = names(v4) = beeCB$Gene
GSEAbias <- lapply(list(v1,v2),GSEAfunc)


ogg11$totBic2=0
ogg11$totBic2[ogg11$abdDE=="conserved worker"]=1
v = ogg11$totBic2

names(v) = ogg11$gene_Mphar
GSEAfunc(v)

factorA$meta = as.factor(apply(factorA[,c(2:6)],1,paste,collapse='_'))
factorB$meta = as.factor(apply(factorB[,c(2:6)],1,paste,collapse='_'))

antCV = averageCV(factorA,antT)
beeCV = averageCV(factorB,beeT)

aCV = data.frame(Gene=names(beeCV[[1]]),CV = beeCV[[1]])
aCVcb = merge(aCV,antCB,by="Gene")
aCVcb2 = merge(aCV,beeRes[[2]],by = "Gene")
aCVcb2$ogg = "no OGG"
aCVcb2$ogg[aCVcb2$ortholog_found] = "paralog present"
aCVcb2$ogg[aCVcb2$OGG_found] = "ortholog present, non-conserved DE"
aCVcb2$ogg[aCVcb2$Gene %in% ogg11$gene_Amel[grepl("conserved",ogg11$abdDE)]] = "conserved DE"
wilcox.test(aCVcb2$CV[aCVcb2$ogg=="conserved DE"],aCVcb2$CV[aCVcb2$ogg=="ortholog present, non-conserved DE"])
ggplot(aCVcb2,aes(x=abdomen,y=CV))+
  geom_boxplot(notch=TRUE)


##Overall network connectivity
beeTn = beeT[rownames(beeT) %in% rownames(bee),]
antTn = antT[rownames(antT) %in% rownames(ant),]
beeTn = log(beeTn + sqrt(1+beeTn^2))
antTn = log(antTn + sqrt(1+antTn^2))

beeCor = cor(t(beeTn))
antCor = cor(t(antTn))
beeCor = data.frame(Gene=rownames(beeCor),kTotal=rowMeans(beeCor),expr=rowMeans(beeTn),Qexpr = rowMeans(beeTn[,grepl("_AQG",colnames(beeTn))]),
                    Wexpr = rowMeans(beeTn[,grepl("_FG|_NG",colnames(beeTn))]))
antCor = data.frame(Gene=rownames(antCor),kTotal=rowMeans(antCor),expr=rowMeans(antTn),Qexpr = rowMeans(antTn[,grepl("_AQG",colnames(antTn))]),
                    Wexpr = rowMeans(antTn[,grepl("_FG|_NG",colnames(antTn))]))
aS = merge(antCor,antCB,by="Gene")
bS = merge(beeCor,beeCB,by="Gene")
aS2 = merge(aS,antRes[[2]],by="Gene")
aS2$Gtype = "non-conserved"
aS2$Gtype[aS2$ortholog_found] = "paralog found"
aS2$Gtype[aS2$OGG_found] = "1-1 OGG"
aS2$Gtype = factor(aS2$Gtype,levels = c("non-conserved","paralog found","1-1 OGG"))
p1 <- ggplot(aS2,aes(x=abdomen,y=Wexpr,fill=Gtype))+
  geom_boxplot()+
  main_theme+
  xlab("abdomen DE")+
  ggtitle("ant")

aS2 = merge(bS,beeRes[[2]],by="Gene")
aS2$Gtype = "non-conserved"
aS2$Gtype[aS2$ortholog_found] = "paralog found"
aS2$Gtype[aS2$OGG_found] = "1-1 OGG"
aS2$Gtype = factor(aS2$Gtype,levels = c("non-conserved","paralog found","1-1 OGG"))

p2 <- ggplot(aS2,aes(x=abdomen,y=Wexpr,fill=Gtype))+
  geom_boxplot(notch=TRUE)+
  main_theme+
  xlab("abdomen DE")+
  ggtitle("bee")

png("~/GitHub/devnetwork/figures/Expr.png",height=2000,width=4000,res=300)
grid.arrange(p1,p2,nrow=1)
dev.off()

ggplot(aS2,aes(x=abdomen,y=kTotal))+
  geom_boxplot(notch=T)

Aps$bic = 0
Aps$bic[Aps$Gene %in% antG] = 1
Bps$bic = 0
Bps$bic[Bps$Gene %in% beeG] = 1

####
#Statistics/numeric results
####
aDE = melt(antRes[[2]][,c(1:6)],id.vars = "Gene")
table(aDE$variable,aDE$value)
bDE = melt(beeRes[[2]][,c(1:6)],id.vars = "Gene")
a=table(bDE$variable,bDE$value)
aDE = melt(antSocRes[[2]][,c(1:4)],id.vars = "Gene")
table(aDE$variable,aDE$value)
bDE = melt(beeSocRes[[2]][,c(1:4)],id.vars = "Gene")
a=table(bDE$variable,bDE$value)
a[,2]+a[,3]

##looking at males
sexDEant = lapply(ant_sexDE,extractBias)
sexDEbee = lapply(bee_sexDE,extractBias)
overQ = c(sum(sexDEant[[3]][[2]] %in% AllPS$Gene.x[AllPS$Gene.y %in% sexDEbee[[3]][[2]]]),length(sexDEant[[3]][[2]]),length(sexDEbee[[3]][[2]]))
overM = c(sum(sexDEant[[3]][[3]] %in% AllPS$Gene.x[AllPS$Gene.y %in% sexDEbee[[3]][[3]]]),length(sexDEant[[3]][[3]]),length(sexDEbee[[3]][[3]]))

Aps$qDE = "non-DE"
Aps$qDE[Aps$Gene %in% sexDEant[[3]][[3]]] = "male"
Aps$qDE[Aps$Gene %in% sexDEant[[3]][[2]]] = "queen"

Bps$qDE = "non-DE"
Bps$qDE[Bps$Gene %in% sexDEbee[[3]][[3]]] = "male"
Bps$qDE[Bps$Gene %in% sexDEbee[[3]][[2]]] = "queen"

conQ = sexDEbee[[3]][[2]][sexDEbee[[3]][[2]] %in% AllPS$Gene.y[AllPS$Gene.x %in% sexDEant[[3]][[2]]]]
conM = sexDEbee[[3]][[3]][sexDEbee[[3]][[3]] %in% AllPS$Gene.y[AllPS$Gene.x %in% sexDEant[[3]][[3]]]]
DmelMale = sexGenes$gene_Amel[sexGenes$FDR < 0.05 & sexGenes$logFC > 0]
DmelF = sexGenes$gene_Amel[sexGenes$FDR < 0.05 & sexGenes$logFC < 0]


p1 <- ggplot(data.frame(x = c(0,1)),aes(x=x))+
  stat_function(fun=function(x) .75*x,color="red",size=2)+
  stat_function(fun=function(x) 1-x,color="green",size=2)+
  stat_function(fun=function(x) 0.5+0.5*x,color="blue",size=2)+
  ylab("trait value")+
  xlab("environment")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20,face="bold"),
        axis.line = element_line(color="black",size=2))
ggsave(p1,file="~/GitHub/devnetwork/figures/dummyPlast.png",height=3,width=3,dpi=300)


f.est <- read.csv("~/Data/Nurse_Larva/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("Gene","f")

prot <- read.csv("~/Data/Nurse_Larva/collectedPAML.csv",sep = "\t",head = F)[,c(1:6)]
colnames(prot) = c('Amel','Sinv','Nvit','A_S','A_N','S_N')

map <- read.table("~/Data/Nurse_Larva/map")
colnames(map) = c("gene_Amel","Amel")
prot = merge(prot,map,by="Amel")
prot$mean = apply(prot[,c(4:6)],1,mean)
prot = merge(prot,ACUogg,by="gene_Amel")

Acb <- merge(prot,antCB,by.x="gene_Mphar",by.y="Gene")
Acb$species = "ant"
Bcb <- merge(prot,beeCB,by.x="gene_Amel",by.y="Gene")
Bcb$species = "bee"


p1C <- ggplot(Acb,aes(x = mean,y=cb_noAbd))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("evolutionary rate")+
  ggtitle("ant")+
  ylab("overall caste bias")+
  #scale_x_log10(breaks = c(0.01,0.1,1))+
  theme(legend.position="none")


p2C <- ggplot(Bcb,aes(x = mean,y=cb_noAbd))+
  geom_hex(bins=70)+
  scale_fill_gradient(low = "blue",high="red")+
  plot2theme+
  ylim(0,2.5)+
  geom_smooth(method="lm",size=1.5,se=FALSE,color="black")+
  xlab("evolutionary rate")+
  ggtitle("bee")+
  ylab("overall caste bias")+
  #scale_x_log10(breaks = c(0.01,0.1,1))+
  theme(legend.position="none")

p <- arrangeGrob(p1C,p2C,ncol=2)
ggsave(p,file="~/GitHub/devnetwork/figures/FigEvol.png",height=6,width=12,dpi=300)

cor.test(Acb$mean,Acb$cb_noAbd,method="spearman")
cor.test(Bcb$mean,Bcb$cb_noAbd,method="spearman")

ogg2 <- read.csv("~/GitHub/devnetwork/data/HymOGG_hym.csv",sep=" ")
t = table(ogg2$OGG)
t = t[t==1]
ogg11 = ogg2[ogg2$OGG %in% names(t),] #get 1-1 orthologs


aEvol = merge(antRes[[2]],f.est,by="Gene")

ggplot(aEvol,aes(x = abdomen,y=f))+
  geom_boxplot(notch=T)

cor.test(aEvol$cb,aEvol$f,method = "spearman")
cor.test(aEvol$cb_Abd,aEvol$f,method = "spearman")
bEvol2 = merge(beeRes[[2]],prot,by.x="Gene",by.y="gene_Amel")
ggplot(bEvol2,aes(x = larva,y=mean))+
  geom_boxplot(notch = T)
bEvol2 = merge(antRes[[2]],prot,by.x="Gene",by.y="gene_Mphar")
ggplot(bEvol2,aes(x = abdomen,y=A_S))+
  geom_boxplot(notch = T)
cor.test(aEvol2$cb_noAbd,aEvol2$A_N,method = "spearman")

bEvol3 = merge(ogg11[!grepl("ant",ogg11$abdDE),],prot,by="gene_Mphar")
ggplot(bEvol3,aes(x=abdDE,y=A_S))+
  geom_boxplot(notch=T)
wilcox.test(bEvol3$mean[bEvol3$abdDE=="conserved queen"],bEvol3$mean[bEvol3$abdDE=="non-conserved/non-DE"])

b = merge(bEvol3,f.est,by.x="gene_Mphar",by.y="Gene")
b = merge(b,beeRes[[2]],by.x="gene_Amel.x",by.y="Gene")
ggplot(b,aes(x=abdomen,y=f))+
  geom_boxplot(notch=T)

DmelSC = merge(sexGenes,ogg11,by="gene_Amel")
key <- read.table("~/GitHub/devnetwork/data/DmelKey.txt") #Generated from Drosophila melanogaster gff file
key2 <- read.table("~/GitHub/devnetwork/data/Dmel_CDStoGene_key.txt")
devel <- read.table("~/Downloads/Devel_IDs.txt")
sex <- read.table("~/Downloads/Sex_IDs.txt")
key3 = merge(key,key2,by="V1")
DmelSC = merge(DmelSC,key3,by.x="Gene",by.y="V2.y")
DmelSC$develSex = "nonDE"
DmelSC$develSex[DmelSC$V2.x %in% devel$V1] = "devel"
DmelSC$develSex[DmelSC$V2.x %in% sex$V1] = "sex"
DmelSC$develSex[DmelSC$V2.x %in% sex$V1 & DmelSC$V2.x %in% devel$V1] = "sex and devel"

DmelSC$develSex = factor(DmelSC$develSex,levels = c("nonDE","sex and devel","devel","sex"))
p <- ggplot(DmelSC[!grepl("ant",DmelSC$abdDE),],aes(x = abdDE,fill=develSex))+
  geom_bar(stat="count",position="fill")+
  ylab("proportion")+
  xlab("abdominal caste bias")
ggsave(p,file = "~/GitHub/devnetwork/figures/SexDevelDmel.png",width=8,height=8)

antDevel2$Gene=rownames(antDevel2)
aD = merge(ogg11,antDevel2,by.x="gene_Mphar",by.y="Gene")

p1 <- ggplot(aD[!grepl("ant",aD$abdDE),],aes(x = abdDE,y=-log(FDR)))+
  geom_boxplot(notch = T)+
  scale_y_log10()+
  ggtitle("devel measured in ant")

beeDevel2$Gene=rownames(beeDevel2)
bD = merge(ogg11,beeDevel2,by.x="gene_Amel",by.y="Gene")

p2 <- ggplot(bD[!grepl("ant",bD$abdDE),],aes(x = abdDE,y=-log(FDR)))+
  geom_boxplot(notch = T)+
  scale_y_log10()+
  ggtitle("devel measured in bee")

p <- arrangeGrob(p1,p2,nrow=1)
ggsave(p,file = "~/GitHub/devnetwork/figures/AbdDevel.png",height=8,width=12,dpi=300)

