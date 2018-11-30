
setwd("~/GitHub/devnetwork/")
load("results/DEtests.RData")
load("results/collectedPhylo.RData")
antConn <- read.csv("results/antConnectivity.csv")
beeConn <- read.csv("results/beeConnectivity.csv")
library(plyr)
library(gridExtra)
library(ggplot2)
library(ppcor)

beeT <- read.table("data/bees.tpm.txt",header=TRUE)
antT <- read.table("data/ants.tpm.txt",header=TRUE)
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

#############
#Table S1
collectDE <- function(d,n1,n2){
  r = as.data.frame(t(apply(d[[2]][,-c(1)],2,function(x){
    c(sum(x!="nonDE"),sum(x==n1),sum(x==n2))
  })))
  colnames(r) = c("total DEGs",paste(n1,"associated"),paste(n2,"associated"))
  r$`stage/tissue` = rownames(r)
  return(r)
}

antC <- collectDE(antRes_allstage,"queen","worker")
beeC <- collectDE(beeRes_allstage,"queen","worker")
antC2 <- collectDE(antRes,"queen","worker")
beeC2 <- collectDE(beeRes,"queen","worker")
antC = rbind(antC,antC2[1,])
beeC = rbind(beeC,beeC2[1,])
antC$`stage/tissue`[rownames(antC) == "larva"] = "larva_overall"
beeC$`stage/tissue`[rownames(beeC) == "larva"] = "larva_overall"
antC$species = "ant"
beeC$species = "honey bee"
C <- rbind(antC,beeC)
C = C[,c(4,5,1:3)]
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(C, rows=NULL, theme=tt)
ggsave(tbl,file = "figures/TableS1.png",height=12,width=10,dpi=300)
write.csv(C,file="figures/TableS1.csv",row.names=FALSE)

############
#Table S2
antC <- collectDE(antSocRes,"nurse","forager")
beeC <- collectDE(beeSocRes,"nurse","forager")
antC$species = "ant"
beeC$species = "honey bee"
C <- rbind(antC,beeC)
C = C[,c(4,5,1:3)]
tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(C, rows=NULL, theme=tt)
ggsave(tbl,file = "figures/TableS2.png",height=12,width=10,dpi=300)
write.csv(C,file="figures/TableS2.csv",row.names=FALSE)

###########
#Table S3
antPlaid <- read.csv("results/antPlaidGenes.csv")
beePlaid <- read.csv("results/beePlaidGenes.csv")
beeAnn <- read.csv("data/Amel_Gene_Names.csv")
antAnn <- read.csv("data/MpharAnn.csv")

antP = merge(antPlaid,antAnn,by="Gene")
antP = antP[,c(1,7,10,11)]
beeP = merge(beePlaid,beeAnn,by="Gene")
beeP = beeP[,c(1,7,10,11)]
antP$abdomen = -antP$abdomen #intially queen was down
beeP$abdomen = -beeP$abdomen #intially queen was down
antP$conn = antP$conn/max(antP$conn)
beeP$conn = beeP$conn/max(beeP$conn)

antP = antP[antP$conn > quantile(antP$conn,0.9)& antP$abdomen > 2,]
antP = antP[order(antP$abdomen,decreasing = T),]
antP = antP[antP$SwissProt!='-',]

beeP = beeP[beeP$conn > quantile(beeP$conn,0.9)& beeP$abdomen > 2,]
beeP = beeP[order(beeP$abdomen,decreasing = T),]
beeP = beeP[beeP$GeneName!='-' & !grepl("uncharacterized",beeP$GeneName),]

colnames(antP) = colnames(beeP) = c("Gene","logFC\nqueen/worker","connectivity","SwissProt")
antP$species = "ant"
beeP$species = "honey bee"
allP = rbind(antP,beeP)

allP$`SwissProt` = gsub("PREDICTED: ","",allP$`SwissProt`)
allP$`logFC\nqueen/worker` = round(allP$`logFC\nqueen/worker`,3)
allP$connectivity = round(allP$connectivity,3)
aP = allP[allP$species=="ant",colnames(allP)!="species"]
bP = allP[allP$species=="honey bee",colnames(allP)!="species"]

aP = aP[!duplicated(aP$SwissProt),]
bP = bP[!duplicated(bP$SwissProt),]

tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(aP[1:30,], rows=NULL, theme=tt)
ggsave(tbl,file = "figures/TableS8.png",height=12,width=14,dpi=300)

tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(bP[1:30,], rows=NULL, theme=tt)
ggsave(tbl,file = "figures/TableS9.png",height=12,width=14,dpi=300)

colnames(allP)[2] = "logFC queen/worker"
write.csv(allP,file="results/TableS3.csv",row.names=F)

########
#Table S4

euclDist <- function(res){
  cb = apply(res[,-c(1)],1,function(x) sqrt(sum(x^2))/length(x))
  r2 = res[,-c(1)]
  cb2 = apply(r2[,!grepl("abdomen",colnames(r2))],1,function(x) sqrt(sum(x^2))/length(x))
  results = data.frame(Gene = res$Gene,cb=cb,cb_noAbd = cb2)
  return(results)
}



partialCor <- function(df,y){
  d = pcor.test(df$cb,df[,colnames(df) == y],df$expr,method="spearman")
  d$variable = y
  d2 = cor.test(df$cb,df[,colnames(df) == y],method="spearman")
  d = cbind(d,data.frame(d2$estimate,d2$p.value))
  d = d[,c(7,8,9,1,2)]
  colnames(d) = c("variable","rho","P-value","rho (partial)","P-value (partial)")
  return(d)
}

antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]])
antSB = euclDist(antSocRes[[1]])
beeSB = euclDist(beeSocRes[[1]])

antCB_expr = euclDist(antRes_allstage[[4]])
beeCB_expr = euclDist(beeRes_allstage[[4]])
antSB_expr = euclDist(antSocRes[[4]])
beeSB_expr = euclDist(beeSocRes[[4]])
antCB = merge(antCB,antCB_expr,by="Gene")
beeCB = merge(beeCB,beeCB_expr,by="Gene")
antSB = merge(antSB,antSB_expr,by="Gene")
beeSB = merge(beeSB,beeSB_expr,by="Gene")

colnames(antCB) = colnames(beeCB) = colnames(antSB) = colnames(beeSB) =c("Gene","cb","cb_noAbd","expr","expr_noAbd")

cor.test(antCB$cb,antCB$expr,method="spearman")
cor.test(antSB$cb,antSB$expr,method="spearman")
cor.test(beeCB$cb,beeCB$expr,method="spearman")
cor.test(beeSB$cb,beeSB$expr,method="spearman")

antCB$type=beeCB$type="caste"
antSB$type=beeSB$type="behavior"
antCB$spec = antSB$spec = "ant"
beeCB$spec = beeSB$spec = "bee"
a = rbind(antCB,antSB)
b = rbind(beeCB,beeSB)

cbAps = merge(a,Aps,by="Gene")
cbBps = merge(b,Bps,by="Gene")

cbAps = merge(antConn,cbAps,by="Gene")
cbBps = merge(beeConn,cbBps,by="Gene")
all = rbind(cbAps,cbBps)


all$kTotal = all$kTotal/max(all$kTotal)
all$psNum = as.numeric(all$psName)
all$ordPs = ordered(all$psName)

lm <- glm(cb ~ kTotal+expr+psName,family="Gamma",data=all[all$cb > 0 & all$type=="caste" & all$spec=="ant",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb ~ kTotal+expr+psName,family="Gamma",data=all[all$cb > 0 & all$type=="caste" & all$spec=="bee",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb ~ kTotal+expr+psName,family="Gamma",data=all[all$cb > 0 & all$type=="behavior" & all$spec=="ant",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb ~ kTotal+expr+psName,family="Gamma",data=all[all$cb > 0 & all$type=="behavior" & all$spec=="bee",])
drop1(lm,.~.,test="Chi")
summary(lm)

##Without abdomen

lm <- glm(cb_noAbd ~ kTotal+expr_noAbd+psName,family="Gamma",data=all[all$cb_noAbd > 0 & all$type=="caste" & all$spec=="ant",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb_noAbd ~ kTotal+expr_noAbd+psName,family="Gamma",data=all[all$cb_noAbd > 0 & all$type=="caste" & all$spec=="bee",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb_noAbd ~ kTotal+expr_noAbd+psName,family="Gamma",data=all[all$cb_noAbd > 0 & all$type=="behavior" & all$spec=="ant",])
drop1(lm,.~.,test="Chi")
summary(lm)

lm <- glm(cb_noAbd ~ kTotal+expr_noAbd+psName,family="Gamma",data=all[all$cb_noAbd > 0 & all$type=="behavior" & all$spec=="bee",])
drop1(lm,.~.,test="Chi")
summary(lm)



levels(cbAps$psName)[1] = levels(cbBps$psName)[1]= "ancient"
cbAps$ordPs = ordered(cbAps$psName)
lm <- glm(cb ~ ordPs,family="Gamma",data=cbAps[cbAps$cb > 0,])
cbBps$ordPs = ordered(cbBps$psName)
lm <- glm(cb ~ ordPs,family="Gamma",data=cbBps[cbBps$cb > 0,])
drop1(lm,.~.,test="Chi")

levels(cbAps$psName)[1] = levels(cbBps$psName)[1]= "ancient"

cor.test(cbAps$kTotal,cbAps$expr,method="spearman")
pcor.test(cbAps$cb,cbAps$kTotal,cbAps$expr,method="spearman")
pcor.test(cbBps$cb,cbBps$kTotal,cbBps$expr,method="spearman")

cbAps$psName = ordered(cbAps$psName)
cbBps$psName = ordered(cbBps$psName)

lm <- glm(log(cb+1) ~ expr + psName + kTotal,data=cbAps,family="gaussian")
drop1(lm,test="Chisq")

lm <- glm(log(cb+1) ~ expr + psName + kTotal,data=cbBps,family="gaussian")
drop1(lm,test="Chisq")

antCB = euclDist(antSocRes[[1]])
beeCB = euclDist(beeSocRes[[1]])

antCB_expr = euclDist(antSocRes[[4]])
beeCB_expr = euclDist(beeSocRes[[4]])
antCB = merge(antCB,antCB_expr,by="Gene")
colnames(antCB) = c("Gene","cb","expr")
beeCB = merge(beeCB,beeCB_expr,by="Gene")
colnames(beeCB) = c("Gene","cb","expr")


cbAps = merge(antCB,Aps,by="Gene")
cbBps = merge(beeCB,Bps,by="Gene")
cbAps = merge(antConn,cbAps,by="Gene")
cbBps = merge(beeConn,cbBps,by="Gene")

cbBps$kTotal = cbBps$kTotal/max(cbBps$kTotal)
cbAps$kTotal = cbAps$kTotal/max(cbAps$kTotal)

levels(cbAps$psName)[1] = levels(cbBps$psName)[1]= "ancient"

cor.test(cbAps$kTotal,cbAps$expr,method="spearman")
pcor.test(cbAps$cb,cbAps$kTotal,cbAps$expr,method="spearman")
pcor.test(cbBps$cb,cbBps$kTotal,cbBps$expr,method="spearman")

cbAps$psName = ordered(cbAps$psName)
cbBps$psName = ordered(cbBps$psName)

lm <- glm(log(cb+1) ~ expr + psName + kTotal,data=cbAps,family="Gamma")
drop1(lm,test="Chisq")

lm <- glm(log(cb+1) ~ expr + psName + kTotal,data=cbBps,family="Gamma")
drop1(lm,test="Chisq")


#########
#GO tables

########
##samples table
df <- data.frame(
  Species = c(rep("ant",28),rep("honey bee",28)),
  Stage = rep(c("egg","L1",rep(c("L2","L3","L4","L5"),each=2),rep("pupa",3),rep("adult",15)),2),
  Tissue = rep(c(rep("whole body",13),rep(c("head","thorax","abdomen"),5)),2),
  Caste = rep(c(rep("N/A",2),rep(c("queen","worker"),4),c("queen","worker","male"),rep(c("queen","worker"),each=6),rep("male",3)),2),
  Type = rep(c(rep("N/A",13),rep("virgin queen",3),rep("mated queen",3),rep("nurse",3),rep("forager",3),rep("N/A",3)),2),
  Replicates = rep(3,56)
)

df$Replicates[c(1,2)] = 6
df$Replicates[39] = 5
df$Replicates[38] = 4

tt <- ttheme_default(colhead=list(fg_params = list(parse=FALSE)),base_size = 16)
tbl <- tableGrob(df, rows=NULL, theme=tt)
ggsave(tbl,file = "figures/SampleTable.png",height=20,width=8,dpi=300)

#########
##Genomes table
df <- read.table("phylostratigraphy/data/chao_codes_edit.txt",head=F,sep="\t")
df = df[,c(1,3)]
colnames(df) = c("Species","NCBI Taxonomy ID")
tbl <- tableGrob(df, rows=NULL, theme=tt)
ggsave(tbl,file = "figures/SpeciesTable.png",height=25,width=8,dpi=300)


