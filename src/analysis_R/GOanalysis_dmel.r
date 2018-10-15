library(topGO)
library(gridExtra)
library(ggplot2)
setwd("~/GitHub/devnetwork/")
load("results/DEtests.RData")
load("results/collectedPhylo.RData")

extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.1 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
}

go <- read.table("data/dmel_ann.txt",sep="\t",header = FALSE,stringsAsFactors = FALSE)

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

GSEAfunc2 <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectConn,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="increasing")
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

antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]])
antSB = euclDist(antSocRes[[1]])
beeSB = euclDist(beeSocRes[[1]])

v1 = antCB$cb
v2 = antSB$cb
names(v1) = names(v2) = antCB$Gene
v3 = beeCB$cb
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
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GOcbTable.png",height=8,width=12,dpi=300)

AcasteRes <- lapply(antTests_oneLarv[c(3:5)],extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],extractBias)

ogg11 = ACUogg
Qc = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[2]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[2]])]
Wc = ogg11$gene_Amel[(ogg11$gene_Amel %in% BcasteRes[[3]][[3]] & ogg11$gene_Mphar %in% AcasteRes[[3]][[3]])]

v1 = v2  = rep(0,nrow(ogg11))
names(v1) = names(v2) = ogg11$gene_Amel
v1[names(v1) %in% Qc] = 1
v2[names(v2) %in% Wc] = 1

spec = c(rep("gene_Amel",2))
cbGO = list(v1,v2)
DmelO <- lapply(seq(1,2),function(i) DmelOrtho(cbGO[[i]],spec[i]))
GO_abd <- lapply(DmelO,GOfunc)


GOcb2 <- lapply(GO_abd,function(x) x[c(1:5),c(1,2,6)])
names = c("conserved queen","conserved worker")
GOcb3 <- lapply(seq(1,2),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3[c(3,4,1,2)])

# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_conDE_Table.png",height=9,width=12,dpi=300)


###Coptying from Figures_IUSSI.R
antPlaid <- read.csv("results/antPlaidGenes.csv")
beePlaid <- read.csv("results/beePlaidGenes.csv")
antG = antPlaid$Gene
beeG = beePlaid$Gene
ogg11=ACUogg
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
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_Bicluster.png",height=9,width=14,dpi=300)


v1 = antPlaid$conn
names(v1) = antPlaid$Gene
v2 = beePlaid$conn
names(v2) = beePlaid$Gene
v1 = DmelOrtho(v1,"gene_Mphar")
v2 = DmelOrtho(v2, "gene_Amel")
ConGO <- lapply(list(v1,v2),GSEAfunc)

GOcb2 <- lapply(ConGO,function(x) x[c(1:5),c(1,2,6)])
names = c("ant","bee")
GOcb3 <- lapply(seq(1,3),function(i) cbind(GOcb2[[i]],test=rep(names[i],5)))
GOcb4 <- do.call(rbind,GOcb3)
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(GOcb4, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_Bicluster2.png",height=9,width=14,dpi=300)

stageGO <- function(df,col,spec){
  c1 = df[,col]
  names(c1) = df$Gene
  DmelO <- DmelOrtho(c1,spec)
  GOres <- GSEAfunc2(DmelO)
  GOcb2 <- GOres[c(1:5),c(1,2,6)]
  GOcb3 <- cbind(GOcb2,test=rep(names(df)[col],5))
  return(GOcb3)
}

nfAnt= ldply(lapply(c(2:4),function(i) stageGO(antSocRes[[3]],i,"gene_Mphar")))
nfBee= ldply(lapply(c(2:4),function(i) stageGO(beeSocRes[[3]],i,"gene_Amel")))
nfGOs = rbind(nfAnt,nfBee)
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(nfAnt, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_nf_ant.png",height=10,width=14,dpi=300)
tbl <- tableGrob(nfBee, rows=NULL, theme=tt)
ggsave(tbl,file = "~/GitHub/devnetwork/figures/GO_nf_bee.png",height=10,width=14,dpi=300)





