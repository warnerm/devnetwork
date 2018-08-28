library(topGO)
library(gridExtra)
library(ggplot2)

extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
}

go <- read.table("~/Downloads/dmel_ann.txt",sep="\t",header = FALSE,stringsAsFactors = FALSE)

new <- list()
for (gene in unique(go$V1)){
  d = go[go$V1 %in% gene,]
  new[[gene]]=as.character(d$V2)
}

selectDE <- function(score){
  return(score == 1)
}

selectConn <- function(score){
  return(score > quantile(score,0.9))
}

DmelOrtho <- function(v,spec){
  v1 = v[names(v) %in% ENDogg[,spec]]
  names(v1) = ENDogg$gene_Dmel[ENDogg[,spec] %in% names(v1)]
  return(v1)
}

#GSEA analysis. Takes in vector of all genes, with a score for social connection strength
GSEAfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectConn,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="decreasing")
  allRes <- GenTable(GOdata,P=resultKS,numChar=100)
  return(allRes) 
}

GOfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectDE,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata,P=resultKS,numChar=100)
  return(allRes) 
}

v1 = antCB$cb_noAbd
v2 = antSB$cb
names(v1) = names(v2) = antCB$Gene
v3 = beeCB$cb_noAbd
v4 = beeSB$cb
names(v3) = names(v4) = beeCB$Gene
spec = c(rep("gene_Mphar",2),rep("gene_Amel",2))
cbGO = list(v1,v2,v3,v4)
DmelO <- lapply(seq(1,4),function(i) DmelOrtho(cbGO[[i]],spec[i]))
GOcb <- lapply(DmelO,GSEAfunc)

GOcb2 <- lapply(GOcb,function(x) x[c(1:5),c(1,2,6)])
names = c("ant caste","ant behavior","bee caste","bee behavior")
GOcb3 <- lapply(seq(1,4),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3)

# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GOcbTable.png",height=12,width=8,dpi=300)

AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

ogg11 = ACUogg
AWBQ = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])] 
AQBW = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])] 
Qc = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])]
Wc = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])]

v1 = v2 = v3 = v4 = rep(0,nrow(ogg11))
names(v1) = names(v2) = names(v3) = names(v4) = ogg11$gene_Amel
v1[names(v1) %in% AWBQ] = 1
v2[names(v2) %in% AQBW] = 1
v3[names(v3) %in% Qc] = 1
v4[names(v4) %in% Wc] = 1

spec = c(rep("gene_Amel",4))
cbGO = list(v1,v2,v3,v4)
DmelO <- lapply(seq(1,4),function(i) DmelOrtho(cbGO[[i]],spec[i]))
GO_abd <- lapply(DmelO,GOfunc)


GOcb2 <- lapply(GO_abd,function(x) x[c(1:5),c(1,2,6)])
names = c("ant queen, bee worker","ant worker, bee queen","conserved queen","conserved worker")
GOcb3 <- lapply(seq(1,4),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3[c(3,4,1,2)])

# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_conDE_Table.png",height=12,width=10,dpi=300)


###Coptying from Figures_IUSSI.R
antQG <- idQG(antPl,9,factorA,antT)
beeQG <- idQG(beePl,7,factorB,beeT)
aG <- commonGenes(antQG,antT)
bG <- commonGenes(beeQG,beeT)
antG = aG$Gene[aG$KeepNum > 0.5]
beeG = bG$Gene[bG$KeepNum > 0.5]
ogg11$antBic = ogg11$beeBic = 0
ogg11$antBic[ogg11$gene_Mphar %in% antG]=1
ogg11$beeBic[ogg11$gene_Amel %in% beeG] = 1
ogg11$totBic = ogg11$antBic + ogg11$beeBic

bicGenesC = bicGenesA = bicGenesB = rep(0,nrow(ogg11))
bicGenesC[ogg11$totBic==2]=1
bicGenesA[ogg11$antBic==1]=1
bicGenesB[ogg11$beeBic==1]=1
names(bicGenesC) = names(bicGenesB) = names(bicGenesA) = ogg11$gene_Amel
bicD = lapply(list(bicGenesA,bicGenesB,bicGenesC),function(x) DmelOrtho(x,"gene_Amel"))
bicGO = lapply(bicD,GOfunc)

GOcb2 <- lapply(bicGO,function(x) x[c(1:5),c(1,2,6)])
names = c("ant","bee","conserved")
GOcb3 <- lapply(seq(1,3),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_Bicluster.png",height=12,width=10,dpi=300)


getConn <- function(tpm,genes,DEtests,DEres){
  tpm = tpm[rownames(tpm) %in% genes,]
  tpm_adj = log(tpm+sqrt(tpm^2+1))
  adj1 = cor(t(tpm_adj))
  conns = data.frame(Gene=rownames(tpm),conn=rowSums(adj1,na.rm = T))
  bDG =DEres[[1]]
  bDG$DEcat = DEres[[2]]$abdomen
  bFDR = data.frame(Gene=rownames(DEtests[["8_gaster"]]),FDR = DEtests[["8_gaster"]]$FDR)
  bDG = merge(bDG,bFDR,by="Gene")
  beeGr = bDG[bDG$Gene %in% genes,]
  beeGrC = merge(beeGr,conns,by="Gene")
  beeGrC$DEcat = factor(beeGrC$DEcat,levels = c("queen","worker","nonDE"))
  return(beeGrC)
}

beeGrC = getConn(beeT,beeG,beeTests,beeRes)
antGrC = getConn(antT,antG,antTests,antRes)

v1 = antGrC$conn
names(v1) = antGrC$Gene
v2 = beeGrC$conn
names(v2) = beeGrC$Gene
v1 = DmelOrtho(v1,"gene_Mphar")
v2 = DmelOrtho(v2, "gene_Amel")
ConGO <- lapply(list(v1,v2),GSEAfunc)

GOcb2 <- lapply(bicGO,function(x) x[c(1:5),c(1,2,6)])
names = c("ant","bee","conserved")
GOcb3 <- lapply(seq(1,3),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_Bicluster.png",height=12,width=10,dpi=300)



