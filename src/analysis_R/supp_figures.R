load("~/GitHub/devnetwork/data/DEtests.RData")
load("~/GitHub/devnetwork/phylo_results/collectedPhylo.RData")
setwd("~/GitHub/devnetwork/src/analysis_R/")
load("../../cluster_results/PlaidResults.RData")

library(RColorBrewer)
SDpalette <- brewer.pal(4,"Set1")

SexPal = c("firebrick2","slateblue4","yellow2","darkgoldenrod2","gray94")

beeT <- read.table("../../data/bees.tpm.txt",header=TRUE)
antT <- read.table("../../data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

#Summarize number of times DE
sumDE <- function(dfDE,type1,type2){
  dfDE$numQueen = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == type1))
  dfDE$numWorker = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == type2))
  d = table(dfDE$numQueen,dfDE$numWorker)
  m = melt(d)
  colnames(m)[c(1,2)] = c(type1,type2)
  return(m)
}

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
    t = cor.test(abs(antD[,i+1]),abs(beeD[,i+1]))
    dAbs[i,2] = t$estimate
    dAbs[i,3] = t$conf.int[1]
    dAbs[i,4] = t$conf.int[2]
  }
  colnames(d) = colnames(dAbs) = c("Stage","cor","c1","c2")
  return(list(d,dAbs))
}


m1 = sumDE(antRes[[2]],"queen","worker")
m2 = sumDE(beeRes[[2]],"queen","worker")

#Add these since there are no genes DE all five times in Apis 
m2E = t(sapply(seq(0,5),function(i) c(queen=i,worker=5,value=0)))
m2Eb = t(sapply(seq(0,5),function(i) c(queen=5,worker=i,value=0)))
m2E = rbind(m2Eb,m2E)
m2 = rbind(m2,m2E)

m1$species = "M. pharaonis"
m2$species = "A. mellifera"
mA = rbind(m1,m2)
mA$species = factor(mA$species,levels= c("M. pharaonis","A. mellifera"))


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
        strip.text = element_text(size=20,face="bold.italic"),
        legend.title = element_text(size = 18,face="bold"),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        panel.border = element_rect(size = 1, color = "black",fill = NA)) 

ggsave(p,file = "~/GitHub/devnetwork/figures/FigS1.png",width = 10,height=8,dpi=300)


m1 = sumDE(antSocRes[[2]],"nurse","forager")
m2 = sumDE(beeSocRes[[2]],"nurse","forager")

m1$species = "M. pharaonis"
m2$species = "A. mellifera"
mA = rbind(m1,m2)
mA$species = factor(mA$species,levels= c("M. pharaonis","A. mellifera"))

#Create heatmap of differential expression (number of times DE for queens and workers)
p <- ggplot(mA,aes(x=nurse,y=forager))+
  geom_tile(aes(fill = value))+
  facet_grid(. ~ species)+
  scale_fill_gradient(name = "number of genes",trans = "log",
                      breaks = c(1,10,100,1000),
                      limits = c(1,10000),
                      labels = c(1,10,100,1000))+
  geom_text(aes(x = nurse,y = forager,label = value),color="white")+
  main_theme+ 
  scale_y_continuous(name = "number of times forager-biased",
                     breaks  = seq(0,5),
                     expand = c(0,0))+
  scale_x_continuous(name = "number of times nurse-biased",
                     breaks = seq(0,5),
                     expand = c(0,0))+
  theme(legend.position = "right",
        axis.line=element_line(color="black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        strip.text = element_text(size=20,face="bold.italic"),
        legend.title = element_text(size = 18,face="bold"),
        strip.background = element_rect(color="black",fill="darkgrey"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        panel.border = element_rect(size = 1, color = "black",fill = NA)) 

ggsave(p,file = "~/GitHub/devnetwork/figures/FigS2.png",width = 10,height=8,dpi=300)

CasteCor <- lfcCor(antRes[[1]],beeRes[[1]])
CasteCor_allStage <- lfcCor(antRes_allstage[[1]],beeRes_allstage[[1]])
BehavCor <- lfcCor(antSocRes[[1]],beeSocRes[[1]])

CasteCor[[1]]$Stage = as.character(CasteCor[[1]]$Stage)
CasteCor[[1]]$Stage[1] = "larva*"
CasteCor_allStage[[1]]$Stage = c("L2","L3","L4","L5","pupa","head","thorax","abdomen")
BehavCor[[1]]$Stage = c("head","thorax","abdomen")
CasteCor[[1]]$Stage = factor(CasteCor[[1]]$Stage,levels = CasteCor[[1]]$Stage)
CasteCor_allStage[[1]]$Stage = factor(CasteCor_allStage[[1]]$Stage,levels = CasteCor_allStage[[1]]$Stage)
BehavCor[[1]]$Stage = factor(BehavCor[[1]]$Stage,levels = BehavCor[[1]]$Stage)


pl <- lapply(list(CasteCor[[1]],CasteCor_allStage[[1]],BehavCor[[1]]), function(x){
  ggplot(x,aes(x = Stage,y=cor))+
    geom_bar(stat="identity")+main_theme+
    xlab("stage/tissue")+
    geom_errorbar(aes(ymin=c1,ymax=c2),width=0.2)+
    ylab("correlation of queen/worker log2 fold-change")+
    theme(plot.margin=unit(c(0.5,1.5,0.5,0.5),"cm"),
          axis.text.x = element_text(hjust=0,angle=-25))
})

pl[[3]] = pl[[3]]+xlab("tissue")+
  ylab("correlation of nurse/forager log2 fold-change")

p <- arrangeGrob(pl[[1]]+annotate("text",x=1,y=0.25,label="A",size=14),
                 pl[[2]]+annotate("text",x=1,y=0.25,label="B",size=14),
                 ncol=2)
ggsave(p,file="~/GitHub/devnetwork/figures/FigS2.png",height=8,width=13.5,dpi=300)
ggsave(pl[[3]],file="~/GitHub/devnetwork/figures/FigS8.png",height=8,width=7,dpi=300)

###Caste vs sex bias
AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)
tissues = c("head","thorax","abdomen")

#Generate plot of log fold change of two DE results, with DE genes highlighted
FCcor <- function(test1,test2){
  FC = merge(test1[[1]],test2[[1]],by = "Gene")
  FC$DE = "nonDE/inconsistent"
  r = cor.test(FC$FC.x,FC$FC.y)
  return(c(r$estimate,c1=r$conf.int[1],c2=r$conf.int[2]))
}

Acor <- ldply(lapply(seq(1,3),function(i){
  FCcor(AcasteRes[[i]],AsexRes[[i]])
}))

Bcor <- ldply(lapply(seq(1,3),function(i){
  FCcor(BcasteRes[[i]],BsexRes[[i]])
}))

Acor$tissue = c("head","thorax","abdomen")
Bcor$tissue = c("head","thorax","abdomen")
Acor$species = "ant"
Bcor$species = "bee"
AllCor = rbind(Acor,Bcor)
AllCor$tissue=factor(AllCor$tissue,levels = c("head","thorax","abdomen"))

p <- ggplot(AllCor,aes(x=tissue,y=cor,fill=species))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(.9),width=0.2)+
  main_theme+
  ylim(0,1)+
  ylab("correlation of caste and sex bias")+
  theme(legend.title = element_text(size=18,face="bold"),
        legend.position = c(0.1,0.9))+
  scale_fill_manual(values = rev(c("goldenrod1","gold4")))

ggsave(p,file = "~/GitHub/devnetwork/figures/FigS3.png",height=6,width=8,dpi=300)

#Sex vs nurse/foraager
AsocRes <- lapply(antSocial,extractBias)
BsocRes <- lapply(beeSocial,extractBias)

Acor <- ldply(lapply(seq(1,3),function(i){
  FCcor(AsocRes[[i]],AsexRes[[i]])
}))

Bcor <- ldply(lapply(seq(1,3),function(i){
  FCcor(BsocRes[[i]],BsexRes[[i]])
}))

Acor$tissue = c("head","thorax","abdomen")
Bcor$tissue = c("head","thorax","abdomen")
Acor$species = "ant"
Bcor$species = "bee"
AllCor = rbind(Acor,Bcor)
AllCor$tissue=factor(AllCor$tissue,levels = c("head","thorax","abdomen"))

p <- ggplot(AllCor,aes(x=tissue,y=cor,fill=species))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(.9),width=0.2)+
  main_theme+
  ylim(-0.2,0.1)+
  ylab("correlation of sex (queen/male)\nand behavior bias (nurse/forager)")+
  theme(legend.title = element_text(size=18,face="bold"),
        legend.position = c(0.1,0.9))+
  scale_fill_manual(values = rev(c("goldenrod1","gold4")))

ggsave(p,file = "~/GitHub/devnetwork/figures/FigS8.png",height=6,width=8,dpi=300)

exprCorPlot <- function(factor,tpm){
  factor$str = as.factor(apply(factor[,c(2,3,5)],1,paste, collapse = "_"))
  f = factor[!duplicated(factor$str),]
  meanExpr <- sapply(fA$str,function(x){
    rowMeans(tpm[,colnames(tpm) %in% factor$sample[factor$str==x]])
  })
  colnames(meanExpr) = f$str
  stageCor <- function(lev,meanExpr,t1,t2){
    ldply(lapply(lev,function(x){
      r = cor.test(
        meanExpr[,grepl(paste(x,t1,sep="_"),colnames(meanExpr))],
        meanExpr[,grepl(paste(x,t2,sep="_"),colnames(meanExpr))]
      )
      c(tissue=x,r$estimate,c1=r$conf.int[1],c2=r$conf.int[2],test=paste(t1,t2,sep="-"))
    }))
  }
  
  lev1 = c("3_larva","4_larva","5_larva","6_larva","7_pupa","8_head","8_mesosoma","8_gaster")
  lev2 = c("7_pupa","8_head","8_mesosoma","8_gaster")
  N1 = c("L2","L3","L4","L5","pupa","head","thorax","abdomen")
  N2 = c("pupa","head","thorax","abdomen")
  QW = stageCor(lev1,meanExpr,"queen","worker")
  QM = stageCor(lev2,meanExpr,"queen","male")
  MW = stageCor(lev2,meanExpr,"worker","male")
  QW$tissue=N1
  MW$tissue=QM$tissue=N2
  allC = do.call(rbind,list(QW,QM,MW))
  allC$tissue = factor(allC$tissue,levels = N1)
  for (i in 2:4){
    allC[,i] = as.numeric(allC[,i])
  }
  allC$test = factor(allC$test,levels = c("queen-worker","queen-male","worker-male"))
  p1 <- ggplot(allC[!grepl("L",allC$tissue),],aes(x = tissue,y=cor,fill=test))+
    geom_bar(stat="identity",position=position_dodge(),alpha=0.8)+
    geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(.9),width=0.2)+
    main_theme+
    ylim(0,1)+
    ylab("expression correlation")+
    theme(legend.title = element_text(size=18,face="bold"),
          legend.position = "left")+
    scale_fill_manual(values = SexPal)
  
  p2 <- ggplot(allC[grepl("L",allC$tissue),],aes(x = tissue,y=cor,fill=test))+
    geom_bar(stat="identity",position=position_dodge(),alpha=0.8)+
    geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(.9),width=0.2)+
    main_theme+
    ylim(0,1)+
    ylab("expression correlation")+
    theme(legend.title = element_text(size=18,face="bold"),
          legend.position = "none")+
    scale_fill_manual(values = SexPal)
  
  return(list(p1,p2))
}

antP <- exprCorPlot(factorA,antT)
beeP <- exprCorPlot(factorB,beeT)

p <- arrangeGrob(
  antP[[1]]+
    theme(plot.margin = unit(c(1,1,0.5,0.5),"cm"),
          legend.text = element_text(size=14)),
  antP[[2]]+
    theme(plot.margin = unit(c(1,0.5,0.5,0.5),"cm")),
  beeP[[1]]+
    theme(plot.margin = unit(c(1,1,0.5,0.5),"cm"),
          legend.text = element_text(size=14)),
  beeP[[2]]+
    theme(plot.margin = unit(c(1,0.5,0.5,0.5),"cm")),
  ncol=2,
  widths = c(0.65,0.35)
)

png("~/GitHub/devnetwork/figures/FigS_expr.png",height=3000,width=3000,res=300)
grid.arrange(p)
grid.text(x=0.2,y=0.98,label="A",gp=gpar(cex=3))
grid.text(x=0.2,y=0.48,label="C",gp=gpar(cex=3))
grid.text(x=0.7,y=0.98,label="B",gp=gpar(cex=3))
grid.text(x=0.7,y=0.48,label="D",gp=gpar(cex=3))
dev.off()



Aps$species = "M. pharaonis"
Bps$species = "A. mellifera"
Aps_FC = merge(Aps,antRes[[1]],by="Gene")
Bps_FC = merge(Bps,beeRes[[1]],by="Gene")
Bps_FC = Bps_FC[!is.na(Bps_FC$psName),]

for (i in 7:11){
  Aps_FC[,i] = Aps_FC[,i] - median(Aps_FC[,i])
  Bps_FC[,i] = Bps_FC[,i] - median(Bps_FC[,i])
}

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

p <- arrangeGrob(pl[[1]]+annotate("text",x=1.2,y=5,label="A. ant",size=14),
                 pl[[2]]+annotate("text",x=1.2,y=5,label="B. bee",size=14),
                 ncol=1)

ggsave(p,file="~/GitHub/devnetwork/figures/FigS5.png",height=12,width=10,dpi=300)

#########

Aps_FC = merge(Aps,antSocRes[[1]],by="Gene")
Bps_FC = merge(Bps,beeSocRes[[1]],by="Gene")
Bps_FC = Bps_FC[!is.na(Bps_FC$psName),]
colnames(Aps_FC)[c(7:9)] = colnames(Bps_FC)[c(7:9)] = c("head","thorax","abdomen")

for (i in 7:9){
  Aps_FC[,i] = Aps_FC[,i] - median(Aps_FC[,i])
  Bps_FC[,i] = Bps_FC[,i] - median(Bps_FC[,i])
}

pl <- lapply(list(Aps_FC,Bps_FC),function(x){
  pm = melt(x,id.vars = c("Gene","ODBgene","OGGacu","ps","psName","species"))
  ggplot(pm,aes(x = variable,y=-value,fill=psName))+
    geom_boxplot(notch=TRUE,outlier.shape = NA)+
    coord_cartesian(ylim = c(-6,6))+
    main_theme+
    scale_fill_manual(values = rev(mypalette),name = "phylostrata")+
    ylab("log2 fold-change\n(nurse/forager)")+
    xlab("tissue")+
    geom_hline(yintercept =0,linetype="dashed")+
    theme(legend.position="right",
          #legend.key.size = unit(0.6,"cm"),
          legend.text = element_text(size=13),
          legend.title = element_text(size=15,face="bold"),
          legend.background = element_blank(),
          plot.margin = unit(rep(1,4),"cm"))
})

p <- arrangeGrob(pl[[1]]+annotate("text",x=1.2,y=5,label="A. ant",size=14),
                 pl[[2]]+annotate("text",x=1.2,y=5,label="B. bee",size=14),
                 ncol=1)

ggsave(p,file="~/GitHub/devnetwork/figures/FigS6.png",height=12,width=10,dpi=300)


###########
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

pl  = list(p1,p2,p3,p4)
pl = lapply(pl,function(x) x + theme(plot.margin = unit(rep(.75,4),"cm")))
lett = c("A","B","C","D")
pl = lapply(seq(1:4),function(i) pl[[i]])
p <- do.call(arrangeGrob,pl)

png("~/GitHub/devnetwork/figures/FigS_bias.png",height=3000,width=3000,res=300)
grid.arrange(p)
grid.text("A",x=0.13,y=0.93,gp=gpar(cex=3))
grid.text("B",x=0.63,y=0.93,gp=gpar(cex=3))
grid.text("C",x=0.13,y=0.43,gp=gpar(cex=3))
grid.text("D",x=0.63,y=0.43,gp=gpar(cex=3))
dev.off()

#Partial correlations with expression
exprA = rowMeans(antT)
exprB = rowMeans(beeT)
exprA = data.frame(Gene = names(exprA),expr = exprA)
exprB = data.frame(Gene = names(exprB),expr = exprB)
antBias = merge(antBias,exprA)
pcor.test(antBias$cb_noAdult.x,antBias$cb.y,antBias$expr)
beeBias = merge(beeBias,exprB)
pcor.test(beeBias$cb_noAdult.x,beeBias$cb.y,beeBias$expr)
aB = merge(antBias,ACUogg,by.x="Gene",by.y="gene_Mphar")
aB = merge(aB,beeBias,by.x="gene_Amel",by.y = "Gene")
pcor.test(aB$cb_noAdult.x.x,aB$cb_noAdult.x.y,aB[,c(15,29)])
pcor.test(aB$cb.y.x,aB$cb.y.y,aB[,c(15,29)])

##########
beeD = rownames(beeDevel2)[beeDevel2$FDR < 0.01]
antD = rownames(antDevel2)[antDevel2$FDR < 0.01]

#overlap of developmental genes
p <- venn.diagram(list("M. pharaonis" = ACUogg$OGGacu[ACUogg$gene_Mphar %in% antD],
                  "A. mellifera" = ACUogg$OGGacu[ACUogg$gene_Amel %in% beeD]),
             filename = NULL,
             imagetype = "png",
             main = "overlap of developmental genes",
             cex=1.5,
             main.cex = 2.2,
             main.pos = c(0.5,0.9),
             cat.cex=1.6,
             height=3000,
             width=5000,
             cat.pos = c(-15,15),
             fontface="bold",
             cat.dist = rep(0.04,2),
             cat.fontface = c("italic","italic"),
             margin = 0.05)


png("~/GitHub/devnetwork/figures/develOverlap.png",height=2000,width=2000,res=300)
grid.draw(p)
dev.off()

a = rbind(c(2279,1589),c(1620,nrow(ACUogg) - 1620 - 2279 - 1589))

develDE <- function(data,devel,conDevel){
  bM = melt(data[[2]],id.vars = "Gene")
  bM$devel = "non-developmental "
  bM$devel[bM$Gene %in% devel]="developmental, one species "
  bM$devel[bM$Gene %in% conDevel]="developmental conserved "
  bM$devel = factor(bM$devel,levels = c("non-developmental ","developmental, one species ","developmental conserved "))
  bM$value = relevel(as.factor(bM$value),ref="nonDE")
  p <- ggplot(bM,aes(x = value,fill=devel))+
    geom_bar(stat="count",position="fill")+
    facet_grid(. ~ variable)+
    main_theme+
    xlab("differential expression")+
    ylab("proportion")+
    scale_fill_manual(values = DE_palette[c(1,2,4)])+
    ylab("proportion of genes")+
    theme(axis.text.x = element_text(angle=-25,hjust=0.1),
          legend.position = "top",
          strip.text = element_text(size=20,face="bold.italic"),
          axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
          legend.title = element_blank(),
          strip.background = element_rect(color="black",fill="darkgrey"),
          plot.margin = margin(.75,.75,.75,.75,"cm"))
  return(p)
}


antConDev = ACUogg$gene_Mphar[ACUogg$gene_Amel %in% beeD & ACUogg$gene_Mphar %in% antD]
beeConDev = ACUogg$gene_Amel[ACUogg$gene_Amel %in% beeD & ACUogg$gene_Mphar %in% antD]

p1 <- develDE(antSocRes,antD,antConDev)
p2 <- develDE(antRes,antD,antConDev)
p3 <- develDE(beeSocRes,beeD,beeConDev)
p4 <- develDE(beeRes,beeD,beeConDev)

grid_arrange_shared_legend <- function(plots,position) {
  g <- ggplotGrob(plots[[1]] + theme(legend.position=position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

png("~/GitHub/devnetwork/figures/DevelProp.png",height=3600,width=6000,res=300)
grid_arrange_shared_legend(list(p1+
                                  theme(legend.text = element_text(size=22)),
                                p2,p3,p4),"bottom")
dev.off()

ogg11 = ACUogg
ogg11$abdDE = "nonDE/inconsistent"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "ant worker, bee queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "ant queen, bee worker"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] = "conserved queen"
ogg11$abdDE[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] = "conserved worker"


ogg11$devel = "non-developmental"
ogg11$devel[ogg11$gene_Mphar %in% antD | ogg11$gene_Amel %in% beeD]="developmental, one species"
ogg11$devel[ogg11$gene_Mphar %in% antD & ogg11$gene_Amel %in% beeD]="developmental conserved"
ogg11$devel = factor(ogg11$devel,levels = c("non-developmental","developmental, one species","developmental conserved"))

ggplot(ogg11,aes(x = abdDE,fill=devel))+
  geom_bar(stat="count",position="fill")

AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)
tissues = c("head","thorax","abdomen")

DEbreakdown <- function(dat,sexDat,devel,t1,t2){
  Ares <- lapply(seq(1,3),function(i){
    d <- data.frame(Gene = sexDat[[1]][[1]]$Gene,DEtype="nonDE ",Caste="nonDE",tissue=tissues[i])
    d$Caste = as.character(d$Caste)
    d$DEtype = as.character(d$DEtype)
    d$Caste[d$Gene %in% dat[[i]][[2]]] = t1
    d$Caste[d$Gene %in% dat[[i]][[3]]] = t2
    d$DEtype[d$Gene %in% devel] = "development "
    d$DEtype[d$Gene %in% c(sexDat[[i]][[2]],sexDat[[i]][[3]])] = "sex "
    d$DEtype[d$Gene %in% devel & d$Gene %in% c(sexDat[[i]][[2]],sexDat[[i]][[3]])] = "sex and development "
    return(d)
  })
  res <- do.call(rbind,Ares)
  res$Caste = relevel(as.factor(res$Caste),ref="nonDE")
  
  res$DEtype = factor(res$DEtype,levels = rev(c("sex and development ","sex ","development ","nonDE ")))
  p <- ggplot(res,aes(x = Caste,fill=DEtype))+
    geom_bar(stat="count",position="fill",alpha=0.8)+
    facet_grid(. ~ tissue)+
    main_theme+
    xlab("differential expression")+
    scale_fill_manual(values = SDpalette)+
    ylab("proportion of genes")+
    theme(axis.text.x = element_text(angle=-25,hjust=0.1),
          legend.position = "top",
          strip.text = element_text(size=20,face="bold.italic"),
          axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
          legend.title = element_blank(),
          strip.background = element_rect(color="black",fill="darkgrey"),
          plot.margin = margin(.75,.75,.75,.75,"cm"))
  return(p)
}

Ac = DEbreakdown(AcasteRes,AsexRes,antD,"queen","worker")
As = DEbreakdown(AsocRes,AsexRes,antD,"nurse","forager")

Bc = DEbreakdown(BcasteRes,BsexRes,beeD,"queen","worker")
Bs = DEbreakdown(BsocRes,BsexRes,beeD,"nurse","forager")

png("~/GitHub/devnetwork/figures/SexAndDevel.png",height=3600,width=6000,res=300)
grid_arrange_shared_legend(list(Ac+
                                  theme(legend.text = element_text(size=22)),
                                As,Bc,Bs),"bottom")
dev.off()

############

antCor = cor(t(antT))^6
beeCor = cor(t(beeT))^6
antConn = data.frame(Gene = rownames(antCor),kTotal = rowSums(antCor))
beeConn = data.frame(Gene = rownames(beeCor),kTotal = rowSums(beeCor))



aM = melt(antRes[[2]],id.vars = "Gene")
bM = melt(beeRes[[2]],id.vars = "Gene")
all = merge(aM,ACUogg,by.x="Gene",by.y="gene_Mphar",all.x=T)
all = merge(all,bM,by.x=c("gene_Amel","variable"),by.y=c("Gene","variable"),all=T)
all$DEspec = "no ortholog"
all$DEspec[all$value.x !="nonDE" & all$value.y != "nonDE"] = "ortholog present, conserved DE"
all$DEspec[all$value.x !="nonDE" & all$value.y == "nonDE" |
             all$value.x =="nonDE" & all$value.y != "nonDE"] = "ortholog present, non-conserved DE"

all$DEspec = factor(all$DEspec,levels = rev(c("ortholog present, conserved DE","ortholog present, non-conserved DE","no ortholog")))
all = merge(all,tau,by.x="gene_Amel",by.y="Gene",all.x=TRUE)
all = merge(all,Aps[,c(1,5)],by = "Gene",all.x=T)
all = merge(all,Bps[,c(1,5)],by.x="gene_Amel",by.y="Gene",all.x=T)
all = merge(all,antConn,by = "Gene",all.x=T)
all = merge(all,beeConn,by.x="gene_Amel",by.y="Gene",all.x=T)
all$beeD = all$antD = "nonDE"
all$beeD[all$gene_Amel %in% beeD] = "DE"
all$antD[all$Gene %in% antD] = "DE"

p3 <- ggplot(all[all$value.y!="nonDE" & !is.na(all$value.y),],aes(x = variable,fill=DEspec,y=tau))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue/stage")+
  theme(legend.position="top")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p1 <- ggplot(all[all$value.x!="nonDE" & !is.na(all$value.x) & !is.na(all$psName.y),],aes(x = variable,fill=DEspec,y=tau))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue/stage")+
  theme(legend.position="top")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p5 <- ggplot(all[all$value.x!="nonDE" & !is.na(all$value.x),],aes(x = variable,fill=DEspec,y=kTotal.x))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue/stage")+
  theme(legend.position="top")+
  scale_y_log10()+
  ylab("log connectivity")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p7 <- ggplot(all[all$value.y!="nonDE" & !is.na(all$value.y),],aes(x = variable,fill=DEspec,y=kTotal.y))+
  geom_boxplot(notch=T)+
  main_theme+
  ylab("log connectivity")+
  xlab("tissue/stage")+
  theme(legend.position="top")+
  scale_y_log10()+
  scale_fill_manual(values = DE_palette[c(1,2,4)])




aM = melt(antSocRes[[2]],id.vars = "Gene")
bM = melt(beeSocRes[[2]],id.vars = "Gene")
all = merge(aM,ACUogg,by.x="Gene",by.y="gene_Mphar",all.x=T)
all = merge(all,bM,by.x=c("gene_Amel","variable"),by.y=c("Gene","variable"),all=T)
all$DEspec = "no ortholog"
all$DEspec[all$value.x !="nonDE" & all$value.y != "nonDE"] = "ortholog present, conserved DE"
all$DEspec[all$value.x !="nonDE" & all$value.y == "nonDE" |
             all$value.x =="nonDE" & all$value.y != "nonDE"] = "ortholog present, non-conserved DE"

all$DEspec = factor(all$DEspec,levels = rev(c("ortholog present, conserved DE","ortholog present, non-conserved DE","no ortholog")))
all = merge(all,tau,by.x="gene_Amel",by.y="Gene",all.x=TRUE)
all = merge(all,antConn,by = "Gene",all.x=T)
all = merge(all,beeConn,by.x="gene_Amel",by.y="Gene",all.x=T)

p4 <- ggplot(all[all$value.y!="nonDE" & !is.na(all$value.y),],aes(x = variable,fill=DEspec,y=tau))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue")+
  theme(legend.position="top")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p2 <- ggplot(all[all$value.x!="nonDE" & !is.na(all$value.x),],aes(x = variable,fill=DEspec,y=tau))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue")+
  theme(legend.position="top")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p6 <- ggplot(all[all$value.x!="nonDE" & !is.na(all$value.x),],aes(x = variable,fill=DEspec,y=kTotal.x))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue")+
  theme(legend.position="top")+
  scale_y_log10()+
  ylab("log connectivity")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])

p8 <- ggplot(all[all$value.y!="nonDE" & !is.na(all$value.y),],aes(x = variable,fill=DEspec,y=kTotal.y))+
  geom_boxplot(notch=T)+
  main_theme+
  xlab("tissue")+
  theme(legend.position="top")+
  scale_y_log10()+
  ylab("log connectivity")+
  scale_fill_manual(values = DE_palette[c(1,2,4)])



p <- arrangeGrob(p1,p2,p3,p4)
p2 <- arrangeGrob(p5,p6,p7,p8)

ggsave(p,file = "~/GitHub/devnetwork/figures/tauDEGs.png",height=10,width=12,dpi=300)
ggsave(p2,file = "~/GitHub/devnetwork/figures/connDEGs.png",height=10,width=12,dpi=300)








