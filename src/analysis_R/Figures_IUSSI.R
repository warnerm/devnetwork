load("~/GitHub/devnetwork/data/DEtests.RData")
load("~/GitHub/devnetwork/phylo_results/collectedPhylo.RData")
library(cowplot)
library(RColorBrewer)
DE_palette = c("#ffffcc","#a1dab4","#41b6c4","#225ea8")
mypalette <- brewer.pal(6,"OrRd")
SexPal = c("firebrick2","slateblue4","yellow2","darkgoldenrod2","gray94")

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
#Calculate tau (tissue specificity) for honey bee genes
calculatetau <- function(factor,expr){
  meanExpr <- lapply(levels(factor$tissue),function(x) 
    rowSums(as.data.frame(expr[,colnames(expr) %in% factor$sample[factor$tissue==x]]))/sum(factor$tissue==x))
  meanExpr <- as.data.frame(do.call(cbind,meanExpr))
  colnames(meanExpr) = levels(factor$tissue)
  geneDeviation = apply(meanExpr,1,function(x){
    sum(ldply(lapply(c(1:12),function(i) 1-x[i]/max(x))))
  })
  
  tau = geneDeviation/(11)
  
  return(list(tau,meanExpr))
}

#Caste-associated genes are also sex-biased
extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
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

p1 <- ggplot(d[d$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  ylim(0,4500)+
  facet_grid(. ~ species)+
  scale_fill_manual(values = DE_palette)+
  ylab("number of caste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.margin = margin(2,2,2,2,"cm"))

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

png("~/GitHub/devnetwork/figures/Fig1a.png",height=2500,width=3500,res=300)
ggdraw()+
  draw_plot(p1+
              theme(legend.text = element_text(size=13),
                    legend.key.width = unit(1,"cm")))+
  draw_plot(p2,x=0.2,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.58,height=0.18,width=0.18)
dev.off()

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
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
beeSocRes[[2]]$ortholog_found = beeSocRes[[2]]$OGG_found = FALSE
beeSocRes[[2]]$ortholog_found[beeSocRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeSocRes[[2]]$OGG_found[beeSocRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeSocRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bMsoc = bM
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

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
  ylab("number of role-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.margin = margin(2,2,2,2,"cm"))

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

png("~/GitHub/devnetwork/figures/Fig1b.png",height=2500,width=3500,res=300)
ggdraw()+
  draw_plot(p1+
              theme(legend.text = element_text(size=13),
                    legend.key.width = unit(1,"cm")))+
  draw_plot(p2,x=0.2,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.58,height=0.18,width=0.18)
dev.off()

##########
##Fig 1d
#########
Aps$species = "M. pharaonis"
Bps$species = "A. mellifera"
Aps_FC = merge(Aps,antRes[[1]],by="Gene")
Bps_FC = merge(Bps,beeRes[[1]],by="Gene")

AllP = rbind(Aps_FC,Bps_FC)
AllP$species = factor(AllP$species,levels=c("M. pharaonis","A. mellifera"))
AllP=droplevels(AllP[!is.na(AllP$psName),])
AllP$psName=as.character(AllP$psName)
AllP$psName[AllP$psName=="bee"|AllP$psName=="ant"]="ant/bee"
AllP$psName=factor(AllP$psName,levels=c("old","insect","hymenoptera","aculeata","ant/bee","novel"))

a = rev(mypalette)

png("~/GitHub/devnetwork/figures/Fig1d.png",width=2500,height=2000,res=300)
ggplot(AllP,aes(x = species,y=-abdomen,fill=psName))+
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
dev.off()

wilcoxLev <- function(vec,levVec,lev1,lev2){
  return(c(lev1,lev2,wilcox.test(vec[levVec==lev1],vec[levVec==lev2])$p.value))
}

res = list()

for (lev1 in levels(Aps_FC$psName)){
  for (lev2 in levels(Aps_FC$psName)){
    res[[lev1]][[lev2]] = wilcoxLev(Aps_FC$abdomen,Aps_FC$psName,lev1,lev2)
  }
}

res = list()

for (lev1 in levels(Bps_FC$psName)){
  for (lev2 in levels(Bps_FC$psName)){
    res[[lev1]][[lev2]] = wilcoxLev(Bps_FC$abdomen,Bps_FC$psName,lev1,lev2)
  }
}

w = sapply(levels(Aps_FC$psName), function(x){
  sapply(levels(Aps_FC$psName), function(y){
    if (x==y) NA
    else wilcoxLev(Aps_FC$abdomen,Aps_FC$psName,x,y)
  })
})

#########
###Figure 2a-c
#########
AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

ogg11 = ACUogg
ogg11$abdDE = "nonDE/inconsistent"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "ant worker, bee queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "ant queen, bee worker"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "conserved queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "conserved worker"

FC = merge(AsexRes[[3]][[1]],AcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Mphar",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "nonDE/inconsistent"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","nonDE/inconsistent"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="nonDE/inconsistent"]=0.8

p1 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal,name = "caste bias")+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  scale_alpha_continuous(guide="none")+
  main_theme+
  ylab("ant caste bias (queen/worker)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="right",
        legend.text = element_text(size=17),
        legend.title = element_text(size=19,face="bold"),
        legend.margin = margin(t=250,r=0,l=0,b=0))

FC = merge(BsexRes[[3]][[1]],BcasteRes[[3]][[1]],by = "Gene")
FC = merge(FC,ogg11,by.x = "Gene",by.y="gene_Amel",all.x=TRUE)
FC$abdDE[is.na(FC$abdDE)] = "nonDE/inconsistent"
FC$abdDE = factor(FC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","nonDE/inconsistent"))
FC$alpha = 0.2
FC$alpha[FC$abdDE!="nonDE/inconsistent"]=0.8

p2 <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  ylab("bee caste bias (queen/worker)")+
  xlab("bee sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="right",legend.title = element_blank())

FCb = merge(BsexRes[[3]][[1]],ogg11,by.x="Gene",by.y="gene_Amel")
FC2 = merge(FCb,AsexRes[[3]][[1]],by.x = "gene_Mphar",by.y= "Gene")
FC2$abdDE = factor(FC2$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","nonDE/inconsistent"))
FC2$alpha = 0.2
FC2$alpha[FC2$abdDE!="nonDE/inconsistent"]=0.8

p3 <- ggplot(FC2,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
  geom_point(aes(fill = abdDE,alpha=alpha),pch=21,color="black",size=2)+
  geom_smooth(method="lm",se=FALSE,color="black")+
  scale_fill_manual(values = SexPal)+
  main_theme+
  ylab("bee sex bias (queen/male)")+
  xlab("ant sex bias (queen/male)")+
  ylim(-10,10)+xlim(-10,10)+
  theme(axis.title.y=element_text(margin = margin(t=0,l=15,r=-5,b=0)))+
  theme(legend.position="right",legend.title = element_blank())

png("~/GitHub/devnetwork/figures/Fig2a_c.png",width=6000,height=2000,res=300)
grid.arrange(grid_arrange_shared_legend(p1,p2,p3,position="right"))
dev.off()

########
##Figure 2d
########
DmelSC = merge(sexGenes,ogg11,by="gene_Amel")
DmelSC$abdDE = factor(DmelSC$abdDE,levels = c("conserved queen","conserved worker","ant worker, bee queen","ant queen, bee worker","nonDE/inconsistent"))
DmelSC$alpha = 0.1
DmelSC$alpha[DmelSC$abdDE!="nonDE/inconsistent"]=0.8

png("~/GitHub/devnetwork/figures/Fig2d.png",width=2000,height=2000,res=300)
ggplot(DmelSC[grepl("conserved",DmelSC$abdDE),],aes(x = abdDE,y=-logFC))+
  #geom_violin(aes(fill=abdDE),alpha=0.8)+
  geom_violin(fill="grey90",trim=FALSE)+
  ylab("logFC (dmel female/dmel male)")+
  geom_jitter(width = 0.1,size=0.5,aes(color=abdDE))+
  geom_boxplot(width=0.05,outlier.shape = NA,fill="black",color="black",notch=TRUE,notchwidth = 0.7)+
  main_theme+
  ylab("fly sex bias (female/male)")+
  ylim(-10,10)+
  xlab("caste bias")+
  scale_fill_manual(values=SexPal,name="caste bias")+
  scale_color_manual(values=SexPal)+
  theme(legend.position="none",
        plot.margin = unit(c(0.5,1.4,0.5,0.5),"cm"))+
  stat_summary(geom = "crossbar", width=0.035, fatten=0, size=0.7,color="white", 
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
  coord_cartesian(ylim = c(-8,8))
dev.off()


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
##Fig 4e
##########
cbAps = merge(antBias2,Aps,by="Gene")
cbBps = merge(beeBias2,Bps,by="Gene")
cbP = rbind(cbAps,cbBps)

cbP$species = factor(cbP$species,levels=c("M. pharaonis","A. mellifera"))
cbP=droplevels(cbP[!is.na(cbP$psName),])
cbP$psName=as.character(cbP$psName)
cbP$psName[cbP$psName=="bee"|cbP$psName=="ant"]="ant/bee"
cbP$psName=factor(cbP$psName,levels=c("old","insect","hymenoptera","aculeata","ant/bee","novel"))
cbP$type = factor(cbP$type,levels = c("caste","behavior"))

p1 <- ggplot(cbP[cbP$type=="caste",],aes(x = species,y=cb_noAdult,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,1.5))+
  scale_fill_manual(values = a,name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall caste bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position="right",
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        legend.margin = margin(t=220,r=0,l=0,b=0), #460 for vertical
        panel.border=element_blank(),
        axis.line.x=element_line("black",size=1),
        axis.line.y=element_line("black",size=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))


p2 <- ggplot(cbP[cbP$type=="behavior",],aes(x = species,y=cb_noAdult,fill=psName))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  main_theme+
  coord_cartesian(ylim=c(0,2.2))+
  scale_fill_manual(values = a,name = "phylostrata",guide=guide_legend(title.position="top",title.hjust = 0.5))+
  ylab("overall behavior bias")+
  theme(axis.text.x = element_text(face="italic"),
        legend.position="right",
        legend.key.size = unit(0.6,"cm"),
        legend.text = element_text(size=15),
        panel.border=element_blank(),
        axis.line.x=element_line("black",size=1),
        axis.line.y=element_line("black",size=1),
        legend.title = element_text(size=17,face="bold"),
        legend.background = element_blank(),
        legend.margin = margin(t=460,r=0,l=-20,b=0),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

png("~/GitHub/devnetwork/figures/Fig4e2.png",width=3500,height=2000,res=300)
grid.arrange(grid_arrange_shared_legend(p1,p2,position="right"))
dev.off()

png("~/GitHub/devnetwork/figures/Fig4e_horiz.png",width=3500,height=1750,res=300)
grid.arrange(p1,p2,ncol=2)
dev.off()

##########
###Fig 5a-d (or supplement)
##########








