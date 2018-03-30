library(edgeR)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)

#Constants
FDR = 0.05

#generalized EdgeR workflow, using quasi-likelihood F-test to identify DEGs
EdgeR_QL <- function(data,design,coef){
  ####Inital steps of standard edgeR analysis
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateDisp(data,design)
  fit <- glmQLFit(dat,design)
  diff <- glmQLFTest(fit, coef=coef) 
  ###Calculate FDR
  out <- topTags(diff,n=Inf,adjust.method="BH")$table
  return(out)
}

#generalized EdgeR workflow, using LR test to identify DEGs
EdgeR <- function(data,design,coef){
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  ##Calculates FDR
  out <- topTags(diff,n=Inf,adjust.method="BH")$table
  return(out)
}

#General function to generate fisher test of a continency table
checkOverlap <- function(dat){
  tbl <- rbind(c(dat[1],dat[2]),
               c(dat[3],dat[4]))
  return(fisher.test(tbl))
}

#filter out based on counts per million reads
filterLowly <- function(counts,cut){ 
  libs = colSums(counts)/1000000
  cpm = counts/libs
  keep = rowSums(cpm > 1) >= ncol(counts)/2
  return(counts[keep,])
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
  return(factors)
}

#Take text file and get data frame we want for analysis
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}

#Put all expression data in one frame by orthologs
getOrthoExpr <- function(bee,ant){
  bee$gene = rownames(bee)
  ant$gene = rownames(ant)
  bee = merge(bee,ogg11,by.x="gene",by.y="gene_Amel")
  ant = merge(ant,ogg11,by.x="gene",by.y="gene_Mphar")
  allD <- merge(bee,ant,by="OGG")
  rownames(allD) = allD$OGG
  allD <- allD[,!grepl("gene",colnames(allD))&!grepl("OGG",colnames(allD))]
  return(allD)
}

#Add tissue_stage bit
editFactor <- function(f){
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male",]) #remove eggs and 1st instars as well as males
  f$tissue_stage = as.factor(do.call(paste,c(f[,c(2,3)],list(sep="_"))))#Want to treat each larval stage separately, as well as each tissue of the adult stage. It's necessary to 
  #do this because the model can't estimate the effect of stage and tissue independently, since there are only adult h/m/g and only larvae of stage 3-6
  relev = levels(f$tissue_stage)
  relev[6:8]=c("8_head","8_mesosoma","8_gaster") #intially alphabetical
  f$tissue_stage = factor(f$tissue_stage,levels=relev)
  return(f)
}

#Treat larva as one stage
editFactor_oneLarv <- function(f){
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male",]) #remove eggs and 1st instars as well as males
  f$tissue_stage = as.factor(do.call(paste,c(f[,c(2,3)],list(sep="_"))))#Want to treat each larval stage separately, as well as each tissue of the adult stage. It's necessary to 
  #do this because the model can't estimate the effect of stage and tissue independently, since there are only adult h/m/g and only larvae of stage 3-6
  relev = levels(f$tissue_stage)
  relev[6:8]=c("8_head","8_mesosoma","8_gaster") #intially alphabetical
  f$tissue_stage = factor(f$tissue_stage,levels=relev)
  f$tissue_stage_oneLarv = f$tissue_stage
  f$tissue_stage[f$tissue == "larva"] = "3_larva"
  f$tissue_stage = droplevels(f$tissue_stage)
  return(f)
}

############
##Loading in, cleaning data
###########
tissues <- c("head","gaster","mesosoma")
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

#load tpm values
beeT <- read.table("~/GitHub/devnetwork/data/bees.tpm.txt",header=TRUE)
antT <- read.table("~/GitHub/devnetwork/data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)

ogg2 <- read.csv("~/GitHub/devnetwork/data/HymOGG_hym.csv",sep=" ")
t = table(ogg2$OGG)
t = t[t==1]
ogg11 = ogg2[ogg2$OGG %in% names(t),] #get 1-1 orthologs

#Generate dataframe of ant and bee expression for each OGG
orthoExpr <- getOrthoExpr(bee,ant)
factorAll <- genFactor(orthoExpr)
factorAll$species="bee"
factorAll$species[grepl("Ant",factorAll$sample)]="ant"

############
###Part 1: Broad signatures of caste - NMDS analysis
############

############
###Part 2: Are there "queen" or "worker" genes across development?
############
#Return list of genes that are DE for caste overall
speciesCasteDE <- function(d,f,factorFun){
  f = factorFun(f)
  d <- d[,colnames(d) %in% f$sample]
  design <- model.matrix(~caste+tissue_stage+colony+caste*tissue_stage,data=f) #Note that we are lumping virgin queens, as well as nurses & foragers
  int <- EdgeR_QL(d,design,(4+length(levels(f$tissue_stage))):(2+length(levels(f$tissue_stage))*2)) 
  genesInt <- rownames(int)[int$FDR < FDR]
  design <- model.matrix(~caste+tissue_stage+colony,data=f) 
  caste <- EdgeR_QL(d,design,2) 
  genesDE <- rownames(caste)[caste$FDR < FDR]
  design <- model.matrix(~caste+tissue_stage+colony,data=f) 
  caste2 <- EdgeR_QL(d[!rownames(d) %in% genesInt,],design,2)  #Remove genes with a significant interaction factor
  genesDE_noInt <- rownames(caste2)[caste2$FDR < FDR]
  return(list(genesDE,genesInt,genesDE_noInt,caste,caste2))
}

#Return list of DE genes at each stage
speciesCasteStage <- function(d,f,factorFun){
  f = factorFun(f)
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue_stage)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue_stage),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste,data=fs) #bee L4 only came from 1 colony
    caste <- EdgeR_QL(ds,design,2)
    results[[lev]]=list(genes=rownames(caste)[caste$FDR<0.05],results=caste)
  }
  return(results)
}

#Take output of 'speciesCasteStage' to identify genes DE at each stage; also genes with consistent logFC
speciesCasteAllStage <- function(stageRes){
  #start with all genes
  genes = geneslfcQ = geneslfcW = rownames(stageRes[[1]][[2]])
  for (i in 1:length(stageRes)){
    genes = genes[genes %in% stageRes[[i]][[1]]]
    geneslfcQ = genes[genes %in% rownames(stageRes[[i]][[2]])[stageRes[[i]][[2]]$logFC < 0]]
    geneslfcW = genes[genes %in% rownames(stageRes[[i]][[2]])[stageRes[[i]][[2]]$logFC > 0]]
  }
  return(list(genes,geneslfcQ,geneslfcW))
}

speciesDEtests <- function(d,f,factorFun){
  OverallDE <- speciesCasteDE(d,f,factorFun)
  Stage <- speciesCasteStage(d,f,factorFun)
  ConsistentDE <- speciesCasteAllStage(Stage)
  return(list(OverallDE,Stage,ConsistentDE))
}

beeTests = speciesDEtests(bee,factorB,editFactor)
beeTests_oneLarv = speciesDEtests(bee,factorB,editFactor_oneLarv)
antTests = speciesDEtests(ant,factorA,editFactor)
antTests_oneLarv = speciesDEtests(ant,factorA,editFactor_oneLarv)

#Finding: There are no "consistent" DE
nGenes_consistentlyDE <- lapply(list(beeTests,beeTests_oneLarv,antTests,antTests_oneLarv),function(x){
  lapply(c(1,2,3),function(j) length(x[[3]][[j]]))
})

#There are, however some genes that are DE overall according to the model. We'll incorporate those below

############
###Part 3: Caste and social toolkits--in each species and the overlap
############
#Return list of DE genes at each stage
speciesSocial <- function(d,f){
  f = droplevels(f[!is.na(f$NF),])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  f$tissue=factor(f$tissue,levels = c("head","mesosoma","gaster"))
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~NF+colony,data=fs)
    res <- EdgeR_QL(ds,design,2)
    results[[lev]]=list(genes=rownames(res)[res$FDR<0.05],results=res)
  }
  return(results)
}

beeSocial <- speciesSocial(bee,factorB)
antSocial <- speciesSocial(ant,factorA)

#Shared social toolkit
sharedGene <- function(b,a){
  b_noOGG = b[!b %in% ogg11$gene_Amel]
  a_noOGG = a[!a %in% ogg11$gene_Mphar]
  aOGG = ogg11$OGG[ogg11$gene_Mphar %in% a]
  bOGG = ogg11$OGG[ogg11$gene_Amel %in% b]
  shared = aOGG[aOGG %in% bOGG]
  
  #Genes with an OGG but not DE in both species
  aOGG_nS = aOGG[!aOGG %in% shared]
  bOGG_nS = bOGG[!bOGG %in% shared]
  
  return(list(shared,a_noOGG,b_noOGG,aOGG_nS,bOGG_nS))
}

social_toolkit <- lapply(seq(1,3), function(i) sharedGene(beeSocial[[i]][[1]],antSocial[[i]][[1]]))

#Caste toolkit across development
caste_toolkit_dev <- lapply(seq(1,8),function(i) sharedGene(beeTests[[2]][[i]][[1]],antTests[[2]][[i]][[1]]))
caste_toolkit_dev_oneLarv <- lapply(seq(1,5),function(i) sharedGene(beeTests_oneLarv[[2]][[i]][[1]],antTests_oneLarv[[2]][[i]][[1]]))

#Caste toolkit, main effect, caste*stage, and main effect with c*s removed
caste_toolkit_overall <- lapply(seq(1,3),function(i) sharedGene(beeTests[[1]][[i]],antTests[[1]][[i]]))
caste_toolkit_overall_oneLarv <- lapply(seq(1,3),function(i) sharedGene(beeTests_oneLarv[[1]][[i]],antTests_oneLarv[[1]][[i]]))

#Correlation of logFC between queens and workers or nurses and foragers at each stage
compLogFC <- function(b,a){
  b$Gene = rownames(b)
  a$Gene = rownames(a)
  b = merge(b,ogg11,by.x="Gene",by.y="gene_Amel")
  a = merge(a,ogg11,by.x="Gene",by.y="gene_Mphar")
  both = merge(a,b,by = "OGG")
  simpCor <- cor.test(both$logFC.x,both$logFC.y)
  absCor <- cor.test(abs(both$logFC.x),abs(both$logFC.y))
  return(list(simpCor,absCor))
}

social_toolkit_lfc <- lapply(seq(1,3), function(i) compLogFC(beeSocial[[i]][[2]],antSocial[[i]][[2]]))

#Caste toolkit across development
caste_toolkit_dev_lfc <- lapply(seq(1,8),function(i) compLogFC(beeTests[[2]][[i]][[2]],antTests[[2]][[i]][[2]]))
caste_toolkit_dev_oneLarv_lfc <- lapply(seq(1,5),function(i) compLogFC(beeTests_oneLarv[[2]][[i]][[2]],antTests_oneLarv[[2]][[i]][[2]]))

#Caste toolkit, main effect and main effect with c*s removed. Note that we can't do caste*stage...
caste_toolkit_overall <- lapply(c(4,5),function(i) compLogFC(beeTests[[1]][[i]],antTests[[1]][[i]]))
caste_toolkit_overall_oneLarv <- lapply(c(4,5),function(i) compLogFC(beeTests_oneLarv[[1]][[i]],antTests_oneLarv[[1]][[i]]))

#Plot number of genes results
plotNgene <- function(res,name){
  numbers <- ldply(lapply(res,function(x) {
    unlist(lapply(x,function(i){
      length(i)
    }))}))
  colnames(numbers) = c("shared","Mphar_noOGG","Amel_noOGG","Mphar_notShared","Amel_notShared")
  numbers$tissue = name
  d = melt(numbers,id.vars = 'tissue')
  p <- ggplot(d,aes(x = tissue, y = value, fill = variable))+
    geom_bar(stat = "identity",position = position_dodge())+
    ylab("Number of DE Genes")
  return(p)
}

#Plot the log fold change comparison results
plotLFC <- function(res,names){
  compiled <- ldply(lapply(unlist(res,recursive = FALSE),function(x) c(x$estimate,x$conf.int,x$p.value)))
  compiled$type = rep(c("signed","absolute value"))
  compiled$name = rep(names,each = 2)
  colnames(compiled)[2:4] = c("ci1","ci2","Pval")
  p <- ggplot(compiled,aes(x = name, y = cor, fill = type))+
    geom_bar(stat = "identity",position = position_dodge())+
    ylab("correlation coefficient")+
    geom_errorbar(aes(ymin = ci1, ymax = ci2),position = position_dodge(width=1),width = 0.5)
  return(p)
}

socPlotLFC <- plotLFC(social_toolkit_lfc,c("head","thorax","abdomen"))
ctdev_lfc_plot <- plotLFC(caste_toolkit_dev_lfc,names(beeTests[[2]]))
ctdev_onelarv_lfc_plot <- plotLFC(caste_toolkit_dev_oneLarv_lfc,names(beeTests_oneLarv[[2]]))
ct_overall_plot <- plotLFC(caste_toolkit_overall,c("Overall","Interaction removed"))
ct_overall_plot_oneLarv <- plotLFC(caste_toolkit_overall_oneLarv,c("Overall","Interaction removed"))

socPlot <- plotNgene(social_toolkit,c("head","thorax","abdomen"))
ctdev_plot <- plotNgene(caste_toolkit_dev,names(beeTests[[2]]))
ctdev_onelarv_plot <- plotNgene(caste_toolkit_dev_oneLarv,names(beeTests_oneLarv[[2]]))
ct_overall_nde_plot <- plotNgene(caste_toolkit_overall,c("Overall","Interaction removed"))
ct_overall_nde_plot_oneLarv <- plotNgene(caste_toolkit_overall_oneLarv,c("Overall","Interaction removed"))

png("~/Writing/Figures/NurseLarva/nDE.png",width=4000,height=5000,res=600)
grid.arrange(socPlot+theme_bw(),
             ctdev_plot+theme_bw(),
             ctdev_onelarv_plot+theme_bw(),ncol = 1)
dev.off()

png("~/Writing/Figures/NurseLarva/logFCfigures.png",width=4000,height=5000,res=600)
grid.arrange(ctdev_lfc_plot+theme_bw(),
             ctdev_onelarv_lfc_plot+theme_bw(),
             ct_overall_plot_oneLarv+theme_bw(),ncol = 1)
dev.off()

png("~/Writing/Figures/NurseLarva/logFCfiguresSoc.png",width=4000,height=5000,res=600)
grid.arrange(socPlotLFC+theme_bw())
dev.off()
############
###Part 4: Caste/social toolkits as components of developmental toolkits
############

############
###Part 5: Phylostrata of caste/social toolkits
############

############
###Part 6: Function of caste/social toolkits
############

#################
##Defining caste toolkits
#################



#Return meta-list of lists of genes that are DE at each stage for a single species
speciesCasteStage <- function(d,f){
  f = editFactor(f)
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue_stage)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue_stage),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste,data=fs) #bee L4 only came from 1 colony
    caste <- EdgeR(ds,design,2)
    results[[lev]]=list(genes=rownames(caste)[caste$FDR<0.05],results=caste)
  }
  return(results)
}

#return list of genes that are DE for caste overall in both species
#Returns caste main effect genes and genes with a caste*stage effect
orthoCaste <- function(d,f){
  f = editFactor(f)
  d <- d[,colnames(d) %in% f$sample]  
  design <- model.matrix(~caste+species+tissue_stage+caste*tissue_stage,data=f) #Note that we are lumping virgin queens, as well as nurses & foragers
  int <- EdgeR(d,design,11:17)
  intGenes <- rownames(int)[int$FDR < FDR]
  design <- model.matrix(~caste+species+tissue_stage,data=f) #Note that we are lumping virgin queens, as well as nurses & foragers
  caste <- EdgeR(d[!rownames(d) %in% intGenes,],design,2) #Remove genes that were found to have an interaction factor
  casteGenes <- rownames(caste)[caste$FDR < FDR]
  f = droplevels(f[factorA$tissue=="larva",])
  design <- model.matrix(~caste+species,data=f)
  casteL <- EdgeR(d[,colnames(d) %in% f$sample],design,2)
  casteLGenes <- rownames(casteL)[casteL$FDR < FDR]
  return(list(main=list(genes=casteGenes,results=caste),int=list(genes=intGenes,results=int),
              casteL = list(genes=casteLGenes,results=casteL)))
}

#Return meta-list of lists of OGGs that are DE at each stage according to a model including species
orthoCasteStage <- function(d,f){
  f = editFactor(f)
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue_stage)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue_stage),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste+species,data=fs) #Don't adjust for colony since since it isn't defined across species
    caste <- EdgeR(ds,design,2)
    results[[lev]]=list(genes=rownames(caste)[caste$FDR<0.05],results=caste)
  }
  return(results)
}
  
beeOverall <- speciesCaste(bee,factorB)
beeStage <- speciesCasteStage(bee,factorB)
antOverall <- speciesCaste(ant,factorA)
antStage <- speciesCasteStage(ant,factorA)
oggOverall <- orthoCaste(orthoExpr,factorAll)
oggStage <- orthoCasteStage(orthoExpr,factorAll)

casteStageOld <- function(d,f){
  f = droplevels(f[f$tissue=="larva"&f$type!="P"&f$stage!=1,])
  results=list()
  for (lev in levels(f$stage)){
    fs = droplevels(f[f$stage==lev,])
    ds <- d[,colnames(d) %in% rownames(fs)]
    if (length(levels(fs$colony)) > 2){
      design <- model.matrix(~caste+colony,data=fs) #Note that we are lumping virgin queens, as well as nurses & foragers
    } else {
      design <- model.matrix(~caste,data=fs) #bee L4 only came from 1 colony
      }
      caste <- EdgeR(ds,design,2)
      results[[lev]]=list(genes=rownames(caste)[caste$FDR<0.05],results=caste)
  }
  return(results)
}
antStageOld <- casteStageOld(antOld,factors)

beeAll <- c(beeStage,beeOverall)
antAll <- c(antStage,antOverall)
oggAll <- c(oggStage,oggOverall)

##################
###Analyze caste toolkit results
##################

#Make bar chart of number of OGG DEs over development
oggSt <- c(oggStage,oggOverall)
res <- matrix(nrow=11,ncol=2)
for (i in 1:length(oggSt)){
  res[i,1]=names(oggSt)[i]
  res[i,2]=length(oggSt[[i]][["genes"]])
}
res = as.data.frame(res)
colnames(res) = c("Comparison","OGG_DE")
levs = c("L2","L3","L4","L5",
         "pupa","head","mesosoma","gaster","Overall Main","CasteXstage")
res$Comparison=levs
res2 = melt(res,id.vars=c("Comparison"))
res2$Comparison=factor(res2$Comparison,levels=levs)
p <- ggplot(res2,aes(x=Comparison,y=as.numeric(value)))+
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=17),
        axis.text.x = element_text(angle=-45,vjust=0.5))+
  ylab("Number of DE genes")

#Make stacked bar charts for the composition of caste DEs across development
#(TRGs, OGG not DE, and OGG DE)
stackedBar <- function(stageDE,overallDE,oggStage,oggOverall,col){
  stageDE <- c(stageDE,overallDE)
  oggSt <- c(oggStage,oggOverall)
  res <- matrix(nrow=10,ncol=4)
  for (i in 1:length(stageDE)){
    res[i,1]=names(stageDE)[i]
    g = stageDE[[i]][["genes"]]
    gO = ogg11$OGG[ogg11[,col] %in% g]
    res[i,2]=sum(oggSt[[i]][["genes"]] %in% gO)
    res[i,3]=length(gO) - as.numeric(res[i,2])
    res[i,4] = length(g)-length(gO)
  }
  res = as.data.frame(res)
  colnames(res) = c("Comparison","OGG_DE","OGG_notDE","noOGG")
  levs = c("L2","L3","L4","L5",
           "pupa","head","mesosoma","gaster","Overall Main","CasteXstage")
  res$Comparison=levs
  res2 = melt(res,id.vars=c("Comparison"))
  res2$Comparison=factor(res2$Comparison,levels=levs)
  p <- ggplot(res2,aes(x=Comparison,y=as.numeric(value),fill=variable))+
    geom_bar(stat="identity")+theme_bw()+
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=17),
          axis.text.x = element_text(angle=-45,vjust=0.5))+
    ylab("Number of DE genes")
  return(list(p,res))
}

antBar <- stackedBar(antStage,antOverall,oggStage,oggOverall,"gene_Mphar")
ggsave(antBar[[1]]+ggtitle("Ant"),file="AntDEs.pdf")
beeBar <- stackedBar(beeStage,beeOverall,oggStage,oggOverall,"gene_Amel")
ggsave(beeBar[[1]]+ggtitle("Bee"),file="BeeDEs.pdf")

#Stacked bar chart of all OGGs found to be DE in the analysis
ogg1A = ogg11[ogg11$gene_Amel %in% rownames(bee)&ogg11$gene_Mphar %in% rownames(ant),] #oggs in analysis
res <- matrix(nrow=10,ncol=4)
fisher <- matrix(nrow=10,ncol=3)
for (i in 1:length(beeAll)){
  res[i,1]=fisher[i,1]=names(beeAll)[i]
  gB = ogg11$OGG[ogg11[,"gene_Amel"] %in% beeAll[[i]][["genes"]]]
  gA = ogg11$OGG[ogg11[,"gene_Mphar"] %in% antAll[[i]][["genes"]]]
  res[i,2]=sum(gB %in% gA)
  res[i,3]=length(gB) - sum(gB %in% gA)
  res[i,4] = length(gA) - sum(gB %in% gA)
  fish=checkOverlap(as.numeric(c(res[i,2:4],nrow(ogg1A))))
  fisher[i,2]=signif(fish$estimate,4)
  fisher[i,3]=signif(fish$p.value,4)
}
res = as.data.frame(res)
fisher=as.data.frame(fisher)
colnames(fisher)=c("Comparison","F","p-value")
colnames(res)=c("Comparison","DEboth","DEapis","DEmphar")
levs = c("L2","L3","L4","L5",
         "pupa","head","mesosoma","gaster","Overall Main","CasteXstage")
res$Comparison=levs
fisher$Comparison=levs
res2 = melt(res,id.vars=c("Comparison"))
res2$Comparison=factor(res2$Comparison,levels=levs)
p <- ggplot(res2,aes(x=Comparison,y=as.numeric(value),fill=variable))+
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=17),
        axis.text.x = element_text(angle=-45,vjust=0.5))+
  ylab("Number of DE genes")
grid.arrange(tableGrob(cbind(res,fisher[,c(2,3)])))

#Take a look at overlapping gaster and interaction effect genes
i=8
gB = ogg11$OGG[ogg11[,"gene_Amel"] %in% beeAll[[i]][["genes"]]]
gA = ogg11$OGG[ogg11[,"gene_Mphar"] %in% antAll[[i]][["genes"]]]
gas = gB[gB %in% gA]
i=10
gB = ogg11$OGG[ogg11[,"gene_Amel"] %in% beeAll[[i]][["genes"]]]
gA = ogg11$OGG[ogg11[,"gene_Mphar"] %in% antAll[[i]][["genes"]]]
int = gB[gB %in% gA]
o = ogg11[ogg11$OGG %in% gas,]

go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
go = go[,c(2:4)]
universe <- as.character(unique(go$gene[go$gene %in% ogg1A$gene_Mphar]))
goFrame=GOFrame(go[go$gene %in% universe,],organism="Monomorium pharaonis")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

GOstat(o$gene_Mphar)

#Make volcano plot for a species DEGs, and color the DEGs that are in the OGG DE list
volcanoToolkit <- function(results1,results2,col1,col2,dT,f){
  upOGG <- rownames(results2)[results2$logFC > 0 & results2$FDR < FDR]
  downOGG <- rownames(results2)[results2$logFC < 0 & results2$FDR < FDR]
  upGene <- ogg11[,col1][ogg11[,col2] %in% upOGG] #Get species names of genes
  downGene <- ogg11[,col1][ogg11[,col2] %in% downOGG]
  f = f[!grepl("L1",f$sample)&!grepl("_E",f$sample)&!grepl("_M",f$sample),]
  dT = dT[rownames(dT) %in% rownames(results1),colnames(dT) %in% f$sample] #Take out 1st instar larvae, eggs, and males
  expr = rowSums(dT)/ncol(dT)
  d = data.frame(expr = expr,gene=names(expr))
  results1$gene=rownames(results1)
  results = merge(results1,d,by="gene")
  results = results[results$FDR < FDR,]
  results$color="black"
  results$color[results$gene %in% upGene]="blue"
  results$color[results$gene %in% downGene]="red"
  p <- ggplot(results,aes(x=log(expr+1),y=logFC))+
        geom_point(color=results$color)+theme_bw()+ylab("logFC queen vs worker")
}

p <- volcanoToolkit(antStage[["8_gaster"]][["results"]],beeStage[["8_gaster"]][["results"]],"gene_Mphar","gene_Amel",antT,factorA[factorA$tissue=="gaster",])
p <- volcanoToolkit(beeStage[["8_gaster"]][["results"]],antStage[["8_gaster"]][["results"]],"gene_Amel","gene_Mphar",beeT,factorB[factorB$tissue=="gaster",])
p <- volcanoToolkit(antOverall[["main"]][["results"]],beeOverall[["main"]][["results"]],"gene_Mphar","gene_Amel",antT,factorA)

#the following one is interesting and weird
p <- volcanoToolkit(beeOverall[["main"]][["results"]],antOverall[["main"]][["results"]],"gene_Amel","gene_Mphar",beeT,factorB)

#Make paired bar chart of expression for the different categories of DEGs
#Could measure expression across castes, or in the caste in which it is higher
pairedBar <- function(DEs,OGGs,dT,f,col){
  f = editFactor(f)
  results <- list()
  for (lev in levels(f$tissue_stage)){
    fs = f[grepl(lev,f$tissue_stage),]
    expr = rowSums(dT[,colnames(dT) %in% fs$sample])/ncol(dT[,colnames(dT) %in% fs$sample])
    expr = data.frame(expr=expr,gene=names(expr))
    res = DEs[[lev]][["results"]]
    res$gene=rownames(res)
    res = merge(res,expr,by="gene")
    res = res[res$FDR < FDR,]
    res$OGG="noOGG"
    res$OGG[res$gene %in% ogg11[,col]]="OGGpresent"
    res$OGG[res$gene %in% ogg11[,col][ogg11$OGG %in% OGGs[[lev]][["genes"]]]]="OGG_DE"
    res$lev = lev
    results[[lev]]=res[,c(7:9)]
  }
  for (lev in c("main","int")){
    expr = rowSums(dT[,colnames(dT) %in% f$sample])/ncol(dT[,colnames(dT) %in% f$sample])
    expr = data.frame(expr=expr,gene=names(expr))
    res = DEs[[lev]][["results"]]
    res$gene=rownames(res)
    res = merge(res,expr,by="gene")
    res = res[res$FDR < FDR,]
    res$OGG="noOGG"
    res$OGG[res$gene %in% ogg11[,col]]="OGGpresent"
    res$OGG[res$gene %in% ogg11[,col][ogg11$OGG %in% OGGs[[lev]][["genes"]]]]="OGG_DE"
    res$lev = lev
    results[[lev]]=res[,c("expr","OGG","lev")] 
  }
  results = ldply(results)
  res = as.data.frame(results)
  colnames(res)= c("Comparison","Expression","OGG")
  p <- ggplot(res,aes(x=Comparison,y=log(Expression+1),fill=OGG))+
    geom_boxplot(notch=TRUE)
  return(list(p,res))
}

p <- pairedBar(antAll,oggAll,antT,factorA,"gene_Mphar")
p <- pairedBar(beeAll,oggAll,beeT,factorB,"gene_Amel")


#Other analyses:
#How old are toolkit genes?
#Are the developmental genes? ask for OGGs as well as within-species
#What functions do they have?
#Could also compare connectivity, or in some way as how central caste toolkit genes are
#Could also be nice to plot logFC of all caste*stage genes
  
#It also would be interesting to plot the expression of MRJP over time

  
######
##Analysis using one model
######
design <- model.matrix(~caste+bigStage+colony+species,data=droplevels(factorAll[factorAll$caste!="male"&factorAll$stage!=1&factorAll$stage!=2,]))
caste <- EdgeR(orthoExpr[,factorAll$caste!="male"&factorAll$stage!=1&factorAll$stage!=2],design,2)

###################
###New section: comparing developmental timecourses
###################


#######
##Heatmaps over time
######
#Filter out 1st and 2nd stage (egg/1st instar) and make a new variable including stage and tissue
filterFactor <- function(f){
  f = f[f$stage!=1&f$stage!=2,] #We'll treat eggs and 1st instar larvae separately since there is no caste
  f$stage_tissue = do.call(paste, c(f[,c("stage","tissue")],list(sep='-')))
  f$stage_tissue = as.factor(f$stage_tissue)
  return(f)
}

#subset data by stage
subData <- function(df,f,lev){
  return(df[colnames(df) %in% f$sample[f$stage_tissue==levels(f$stage_tissue)[lev]]])
}

#calculate correlation matrix across stages for two separate dataframes
corStage <- function(df1,df2,f1,f2){
  f1 = filterFactor(f1)
  f2 = filterFactor(f2)
  corMat <- matrix(nrow=8,ncol=8)
  for (i in 1:8){
    for (j in 1:8){
      d1 = subData(df1,f1,i)
      d2 = subData(df2,f2,j)
      corMat[i,j] = cor((rowSums(d1)/ncol(d1)),(rowSums(d2)/ncol(d2)))
    }
  }
  return(corMat)
}

makePlot <- function(corMat,labx,laby){
  m = as.data.frame(corMat)
  colnames(m) = rownames(m) = c("L2","L3","L4","L5","pupa","gaster","head","mesosoma")
  m$StageX = factor(rownames(m),levels=colnames(m))
  m1 = melt(m,id.vars="StageX")
  p <- ggplot(m1,aes(variable,StageX))+
    geom_tile(aes(fill=value))+
    xlab(labx)+ylab(laby)+
    scale_fill_continuous(na.value="white",name="r",guide="colorbar")+
    guides(fill=guide_colorbar(nbin=100))+
    scale_y_discrete(limits=rev(levels(m1$StageX)),position="right")+
    scale_x_discrete(position="top")+
    theme(legend.position="bottom",
          legend.key.width=unit(1,"cm"),
          panel.background = element_rect(fill="white"),
          axis.text=element_text(size=11),
          axis.title=element_text(size=17),
          axis.title.x=element_text(margin=margin(t=20,r=0,b=20,l=0)),
          axis.title.y=element_text(margin=margin(t=0,r=0,b=0,l=20)),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm"))
  return(p)
}

antbee = corStage(orthoExpr[,grepl("Bee",colnames(orthoExpr))],orthoExpr[,grepl("Ant",colnames(orthoExpr))],
                  factorAll[grepl("Bee",rownames(factorAll)),],factorAll[grepl("Ant",rownames(factorAll)),])

antant = corStage(orthoExpr[,grepl("Ant",colnames(orthoExpr))],orthoExpr[,grepl("Ant",colnames(orthoExpr))],
                  factorAll[grepl("Ant",rownames(factorAll)),],factorAll[grepl("Ant",rownames(factorAll)),])


makePlot(beeQW,"Bee","Ant")

beeQW = corStage(beeT[,factorB$caste=='queen'],beeT[,factorB$caste=='worker'],factorB[factorB$caste=='queen',],factorB[factorB$caste=='worker',])
antQW = corStage(antT[,factorA$caste=='queen'],antT[,factorA$caste=='worker'],factorA[factorA$caste=='queen',],factorA[factorA$caste=='worker',])


#########Functions for the analysis
#############
###Devel toolkit
casteDev <- function(fdr,caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(rownames(out)[out$FDR < fdr])
}

genDevTool <- function(fdr,factors,data){
  fQueen <- factors[factors$stage==1|factors$stage==2,]
  fQueen$caste="queen"
  fQueen$sample=paste(fQueen$sample,"_QueenCopy",sep="")
  rownames(fQueen)=fQueen$sample
  factors <- rbind(factors,fQueen) #Add copies of egg and L1 samples for dev toolkit definition
  dataQ <- data[,grepl("L1|_E",colnames(data))]
  colnames(dataQ)=paste(colnames(dataQ),"_QueenCopy",sep="")
  data <- cbind(data,dataQ)#Add copies of egg and L1 samples for dev toolkit definition
  QGenes <- casteDev(fdr,"queen",data,factors)
  WGenes <- casteDev(fdr,"worker",data,factors)
  return(QGenes[QGenes %in% WGenes]) #return genes that are differentially expressed across stages in queens and workers
}

##Entire workflow with different fdr values
workflow <- function(fdr){
  BeeDev <- genDevTool(fdr,factorB,bee)
  AntDev <- genDevTool(fdr,factorA,ant)
  
  BeeDevOGG <- ogg2$OGG[ogg2$gene_Amel %in% BeeDev]
  AntDevOGG <- ogg2$OGG[ogg2$gene_Mphar %in% AntDev]
  TwoSpecDev = BeeDevOGG[BeeDevOGG %in% AntDevOGG]
  
  DevList <- list(BeeDev,AntDev,TwoSpecDev)
  nGene <- list(nrow(bee),nrow(ant),length(commonOGG))
  
  beeCaste <- list()
  antCaste <- list()
  TwoSpecCaste <- list()
  
  for (tissue in tissues){
    beeCaste[[tissue]]=tissueCaste(factorB,bee,tissue,fdr)
    antCaste[[tissue]]=tissueCaste(factorA,ant,tissue,fdr)
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% beeCaste[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% antCaste[[tissue]]]
    TwoSpecCaste[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  Castelist <- list(beeCaste,antCaste,TwoSpecCaste)
  
  beeNF <- list()
  antNF <- list()
  TwoSpecNF <- list()
  
  for (tissue in tissues){
    beeNF[[tissue]]=tissueSocial(factorB,bee,tissue,fdr)
    antNF[[tissue]]=tissueSocial(factorA,ant,tissue,fdr)
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% beeNF[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% antNF[[tissue]]]
    TwoSpecNF[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  NFlist <- list(beeNF,antNF,TwoSpecNF)
  
  #######
  ##Hypothesis tests
  ######Make table with # of genes in each of the lists above (not developmental), plus number of genes in analysis (in parentheses)
  ######Then add number of genes in each overlap comparison
  ######Rows are head, mesosoma, and gaster for both caste and social
  ######Can do a quick chi square, but also do bootstrapping for confidence intervals based on scrambled data
  test_list <- list(Castelist,NFlist)
  
  table <- matrix(nrow = 8, ncol = 9)
  for (i in 1:2){
    for (j in 1:3){
      for (k in 1:3){
        table[(i-1)*3+j,k] = length(test_list[[i]][[k]][[tissues[j]]])
      }
    }
  }
  
  for (i in 1:3){
    table[7,i] = length(DevList[[i]])
    table[7,i+3] = length(DevList[[i]])
    table[7,i+6] = length(DevList[[i]])
    table[8,i] = nGene[[i]]
    table[8,i+3] = nGene[[i]]
    table[8,i+6] = nGene[[i]]
  }
  
  
  
  for (i in 1:3){
    for (j in 1:2){
      for (k in 1:3){
        table[(j-1)*3+k,i+3] = sum(DevList[[i]] %in% test_list[[j]][[i]][[tissues[k]]])
        chi = rbind(c(table[(j-1)*3+k,i+3],
                      table[(j-1)*3+k,i]-table[(j-1)*3+k,i+3]),
                    c(table[7,i] - table[(j-1)*3+k,i+3],
                      nGene[[i]] - table[7,i] - table[(j-1)*3+k,i] + table[(j-1)*3+k,i+3]))
        table[(j-1)*3+k,i+6] = signif(fisher.test(chi,alternative="greater")$p.value,digits=3)
      }
    }
  }
  
  rownames(table) <- c(apply(expand.grid(tissues,c("Caste","NF")), 1, paste, collapse="."),"DevelToolkit","nGene_Analysis")
  colnames(table) <- c("Amel","Mpharaonis","TwoSpecies","Amel_Devel",
                       "Mpharaonis_Devel","TwoSpecies_Devel",
                       "Amel_fisherP","Mpharaonis_fisherP","TwoSpecies_fisherP")
  
  png(paste("~/GitHub/devnetwork/results/OverlapTable",fdr,"fdr",".png",sep=""),width=4000,height=1000,res=300)
  grid.table(table)
  dev.off()
  
  write.table(table,file=paste("~/GitHub/devnetwork/results/OverlapTable",fdr,".txt",sep=""),row.names=TRUE)
  return(list(test_list,DevList))
}


##########
##Caste toolkit
tissueCaste <- function(factors,data,tissue,fdr,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste!="male"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$caste = sample(f$caste,length(f$caste),replace=F)
  }
  design <- model.matrix(~caste+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < fdr])
}

#############
##Social toolkit
tissueSocial <- function(factors,data,tissue,fdr,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste=="worker"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$NF = sample(f$NF,length(f$NF),replace=F)
  }
  design <- model.matrix(~NF+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(rownames(out)[out$FDR < fdr])
}

#constants
tissues <- c("head","gaster","mesosoma")

BeeOGG <- ogg2$OGG[ogg2$gene_Amel %in% rownames(bee)]
AntOGG <- ogg2$OGG[ogg2$gene_Mphar %in% rownames(ant)]
commonOGG <- BeeOGG[BeeOGG %in% AntOGG]

#############
##Running program
#############
fdr_0.05 <- workflow(0.05)
fdr_0.1 <- workflow(0.1)
fdr_0.3 <- workflow(0.3)

############
##Loading bootstrapped results
###########
casteBoot <- read.csv("~/GitHub/devnetwork/data/data.processed/casteBoot.csv")
socialBoot <- read.csv("~/GitHub/devnetwork/data/data.processed/socialBoot.csv")

dfBoot <- function(test,fdr,bootRes,trueRes){
  bootRes <- bootRes[-c(1),-c(1)]
  res <- data.frame(type=rep("bee",9),trueDE=rep(length(trueRes[[1]][["head"]]),9),
                    bootDE=NA,
                    bootC95=NA,
                    bootC99=NA)
  res$type = colnames(bootRes)
  for (i in 1:9){
    res$bootDE[i] = round(mean(bootRes[,i],na.rm=TRUE),3)
    res$bootC95[i] = round(quantile(bootRes[,i],0.95,na.rm=TRUE),3)
    res$bootC99[i] = round(quantile(bootRes[,i],0.99,na.rm=TRUE),3)
  }
  for (i in 1:3){
    for (j in 1:3){
      res$trueDE[(i-1)*3+j] = length(trueRes[[i]][[tissues[[j]]]])
    }
  }
  return(res)
}

caste0.05 <- dfBoot("caste",0.05,casteBoot,fdr_0.05[[1]][[1]])
rownames(caste0.05)=caste0.05$type
plotDf <- melt(caste0.05,id.vars=c("type","bootC95","bootC99"))
plotDf$bootC95[plotDf$variable=="trueDE"]=NA
plotDf$bootC99[plotDf$variable=="trueDE"]=NA
ggplot(plotDf,aes(x=type,y=value,fill=variable))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymax=bootC99,ymin=value),position=position_dodge())

caste0.05=caste0.05[,-c(1)]
png("~/GitHub/devnetwork/results/OverlapTableBootCaste.png",width=2000,height=1000,res=300)
grid.table(caste0.05)
dev.off()


caste0.05 <- dfBoot("caste",0.05,socialBoot,fdr_0.05[[1]][[2]])
rownames(caste0.05)=caste0.05$type
caste0.05=caste0.05[,-c(1)]
png("~/GitHub/devnetwork/results/OverlapTableBootSocial.png",width=2000,height=1000,res=300)
grid.table(caste0.05)
dev.off()

###########
##Trying analysis by looking at correlation of F values between Dev toolkit and caste/social toolkit
##########
EdgeR <- function(data,design,coef){
  ####Inital steps of standard edgeR analysis
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  ###Calculate FDR
  out <- topTags(diff,n=Inf,adjust.method="BH")$table
  return(out)
}


casteDev <- function(caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(out)
}

genDevTool <- function(factors,data){
  fQueen <- factors[factors$stage==1|factors$stage==2,]
  fQueen$caste="queen"
  fQueen$sample=paste(fQueen$sample,"_QueenCopy",sep="")
  rownames(fQueen)=fQueen$sample
  factors <- rbind(factors,fQueen) #Add copies of egg and L1 samples for dev toolkit definition
  dataQ <- data[,grepl("L1|_E",colnames(data))]
  colnames(dataQ)=paste(colnames(dataQ),"_QueenCopy",sep="")
  data <- cbind(data,dataQ)#Add copies of egg and L1 samples for dev toolkit definition
  worker <- casteDev("worker",data,factors)
  queen <- casteDev("queen",data,factors)
  return(list(worker,queen))
}

getOGG <- function(out,merge.col){
  out$gene = rownames(out)
  out <- merge(out,sharedOGG,by.x="gene",by.y=merge.col)
  return(out)
} 

sharedOGG <- ogg2[ogg2$gene_Mphar %in% rownames(ant) & ogg2$gene_Amel %in% rownames(bee),]
#remove duplicates
sharedOGG$num=1
for (i in 1:nrow(sharedOGG)){
  num = nrow(sharedOGG[sharedOGG$OGG==sharedOGG$OGG[i],])
  sharedOGG$num[i] = num
}
sharedOGG <- sharedOGG[sharedOGG$num==1,]

beeDev <- genDevTool(factorB,bee)
antDev <- genDevTool(factorA,ant)
beeDev_OGG <- lapply(beeDev,getOGG,merge.col="gene_Amel")
antDev_OGG <- lapply(antDev,getOGG,merge.col="gene_Mphar")
Dev_OGG <- c(beeDev_OGG,antDev_OGG)

venn = matrix(nrow=nrow(beeDev_OGG[[1]]),ncol=4)
Fmat = matrix(nrow=nrow(beeDev_OGG[[1]]),ncol=4)
for (i in 1:4){
  for (j in 1:nrow(sharedOGG)){
    row = Dev_OGG[[i]][Dev_OGG[[i]]$OGG %in% sharedOGG$OGG[j],]
    Fmat[j,i] = row$LR
    if (row$FDR < FDR){
      venn[j,i] = 1
    } else {
      venn[j,i] = 0
    }
  }
}

x <- vennCounts(venn)
png("~/GitHub/devnetwork/results/devVenn.png")
vennDiagram(x,names=c("beeW","beeQ","antW","antQ"),main="Overlap of developmental DE genes")
dev.off()

corMat = matrix(nrow=4,ncol=4)
colnames(corMat) = rownames(corMat) = c("beeW","beeQ","antW","antQ")
for (i in 1:4){
  for (j in 1:4){
    corMat[i,j] = round(cor(Fmat[,i],Fmat[,j]),3)
  }
}
mytheme <- gridExtra::ttheme_default(
  core = list(padding=unit(c(1, 1), "mm"))
)
t <- tableGrob(corMat,theme=mytheme)
png("~/GitHub/devnetwork/results/devCor.png")
grid.arrange(t,top=textGrob("pearson correlation of LR for each gene
            where LR is testing the affect of stage on expression.
            (i.e. whether or not a gene changes expression across pre-pupal development)
            P <<< 0.001 in all cases",vjust=3))
dev.off()

bVenn <- matrix(data=0,nrow=nrow(venn),ncol=2)
bVenn[,1][venn[,1]+venn[,2]==2]=1
bVenn[,2][venn[,3]+venn[,4]==2]=1
x <- vennCounts(bVenn)

chisq.test(rbind(c(x[3,3],x[2,3]),c(x[1,3],x[0,3])))
png("~/GitHub/devnetwork/results/devVenn2.png")
vennDiagram(x,names=c("Amel Dev","Mphar Dev"))
dev.off()

##########
##Caste toolkit
tissueCaste <- function(factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste!="male"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$caste = sample(f$caste,length(f$caste),replace=F)
  }
  design <- model.matrix(~caste+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(out)
}

#############
##Social toolkit
tissueSocial <- function(factors,data,tissue,scramble=FALSE){
  counts <- data[,factors$stage==8&factors$tissue==tissue&factors$caste=="worker"]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  if (scramble){ #mix up caste labels
    f$NF = sample(f$NF,length(f$NF),replace=F)
  }
  design <- model.matrix(~NF+colony,data=f)
  out <- EdgeR(counts,design,2)
  return(out)
}

casteList <- list(list())
socList <- list(list())
for (tissue in tissues){
  casteList[["bee"]][[tissue]]=tissueCaste(factorB,bee,tissue)
  socList[["bee"]][[tissue]]=tissueSocial(factorB,bee,tissue)
  casteList[["ant"]][[tissue]]=tissueCaste(factorA,ant,tissue)
  socList[["ant"]][[tissue]]=tissueSocial(factorA,ant,tissue)
}

numDEframe <- function(list){
  beeOGG <- lapply(list[["bee"]],getOGG,merge.col="gene_Amel")
  antOGG <- lapply(list[["ant"]],getOGG,merge.col="gene_Mphar")
  fr = matrix(nrow=3,ncol=3)
  corMat = matrix(nrow=3,ncol=4)
  for (i in 1:3){
    over = merge(beeOGG[[tissues[i]]],antOGG[[tissues[i]]],by="OGG")
    fr[i,1]=sum(beeOGG[[tissues[i]]]$FDR < FDR)
    fr[i,2]=sum(antOGG[[tissues[i]]]$FDR < FDR)
    fr[i,3]=sum(over$FDR.x < 0.05 & over$FDR.y < 0.05)
    mat = cor.test(over$logFC.x,over$logFC.y)
    corMat[i,1] = signif(mat$estimate,4)
    corMat[i,2] = signif(mat$p.value,4)
    mat = cor.test(abs(over$logFC.x),abs(over$logFC.y))
    corMat[i,3] = signif(mat$estimate,4)
    corMat[i,4] = signif(mat$p.value,4)
  }
  rownames(fr)= rownames(corMat) = tissues
  colnames(fr) = c("bee","ant","overlap")
  colnames(corMat) = c("R","pVal","R_abs","pVal_abs")
  return(list(fr,corMat))
}

cDE <- numDEframe(casteList)
sDE <- numDEframe(socList)

res = c(cDE,sDE)
l = lapply(res,tableGrob)
png("~/GitHub/devnetwork/results/casteTables.png",height=1000,width=3000,res=300)
grid.arrange(l[[1]],l[[3]],l[[2]],l[[4]])
dev.off()

devGenes <- function(dev){
  dev2 <- merge(dev[[1]],dev[[2]],by="OGG")
  return(dev2$OGG[dev2$FDR.x < 0.05 & dev2$FDR.y < 0.05])
}

devLR <- function(dev){
  dev2 <- merge(dev[[1]],dev[[2]],by="OGG")
  dev2$meanLR = (dev2$LR.x+dev2$LR.y)/2
  return(dev2)
}
dev = list(beeDev_OGG,antDev_OGG)
devGene <- lapply(dev,devGenes)
devLRlist <- lapply(dev,devLR)

numDEoverlap <- function(list,test){
  beeOGG <- lapply(list[["bee"]],getOGG,merge.col="gene_Amel")
  antOGG <- lapply(list[["ant"]],getOGG,merge.col="gene_Mphar")
  nDE <- matrix(nrow=3,ncol=9)
  corMat <- matrix(nrow=3,ncol=4)
  for (i in 1:3){
    beeG = beeOGG[[tissues[i]]][beeOGG[[tissues[i]]]$FDR < FDR,]
    antG = antOGG[[tissues[i]]][antOGG[[tissues[i]]]$FDR < FDR,]
    overG = beeG$OGG[beeG$OGG %in% antG$OGG]
    nDE[i,1] = length(devGene[[1]])
    nDE[i,2] = nrow(beeG)
    nDE[i,3] = sum(beeG$OGG %in% devGene[[1]])
    nDE[i,4] = length(devGene[[2]])
    nDE[i,5] = nrow(antG)
    nDE[i,6] = sum(antG$OGG %in% devGene[[2]])
    nDE[i,7] = sum(devGene[[1]] %in% devGene[[2]])
    nDE[i,8] = length(overG)
    nDE[i,9] = sum(overG %in% devGene[[1]][devGene[[1]] %in% devGene[[2]]])
    mat = cor.test(devLRlist[[1]]$meanLR,beeOGG[[tissues[i]]]$LR)
    corMat[i,1] = signif(mat$estimate,4)
    corMat[i,2] = signif(mat$p.value,4)
    mat = cor.test(devLRlist[[2]]$meanLR,antOGG[[tissues[i]]]$LR)
    corMat[i,3] = signif(mat$estimate,4)
    corMat[i,4] = signif(mat$p.value,4)
  }
  rownames(nDE) = rownames(corMat) = tissues
  colnames(nDE) =  apply(expand.grid(c("devel",test,paste("devel_",test,sep="")),c("bee","ant","overlap")), 1, paste, collapse=".")
  colnames(corMat) = c("bee_R","bee_pval","ant_R","ant_pval")
  return(list(nDE,corMat))
}

cOver <- numDEoverlap(casteList,"caste")
sOver <- numDEoverlap(socList,"social")

tbl = list(cOver[[1]],sOver[[1]])
tbls <- lapply(tbl,tableGrob)

png("~/GitHub/devnetwork/results/overTables.png",height=1500,width=4000,res=300)
grid.arrange(tbls[[1]],tbls[[2]])
dev.off()

##############
##Comparing to genes defined as developmental toolkit in D. melanogaster
#############
dmel <- read.table("~/GitHub/devnetwork/data/DmelDevelOGG.txt")
oggDev = droplevels(ogg3[ogg3$OGG %in% dmel$V1,])
t = table(oggDev$OGG)
t = t[t==1]
oggDevF = oggDev[oggDev$OGG %in% names(t),] #Filter for 1-1-1 orthologs (developmental genes)

t = table(ogg3$OGG)
t = t[t==1]
ogg111 = ogg3[ogg3$OGG %in% names(t),] #Filter for 1-1-1 orthologs (all genes)
ogg111 = ogg111[ogg111$gene_Mphar %in% rownames(antDev[[1]]) &
                  ogg111$gene_Amel %in% rownames(beeDev[[1]]),] #Filter for oggs in the analysis


beeDW = rownames(beeDev[[1]])[beeDev[[1]]$FDR < FDR]
antDW = rownames(antDev[[1]])[antDev[[1]]$FDR < FDR]
beeDQ = rownames(beeDev[[2]])[beeDev[[2]]$FDR < FDR]
antDQ = rownames(antDev[[2]])[antDev[[2]]$FDR < FDR]

beeD = beeDW[beeDW %in% beeDQ]
antD = antDW[antDW %in% antDQ]

beeD_ogg = ogg111$OGG[ogg111$gene_Amel %in% beeD]
antD_ogg = ogg111$OGG[ogg111$gene_Mphar %in% antD]
DevOGG = beeD_ogg[beeD_ogg %in% antD_ogg]

MasterDev = DevOGG[DevOGG %in% oggDevF$OGG] #OGGs identified as dev toolkit by transcriptomic and Dmel annotation

a = rbind(c(length(MasterDev),nrow(oggDevF) - length(MasterDev)),
          c(length(DevOGG) - length(MasterDev),
            nrow(ogg111) - nrow(oggDevF) - length(DevOGG) + length(MasterDev)))

fisher.test(a,alternative="greater") #They is less association than expected
oggDevF[oggDevF$OGG %in% MasterDev,]

##############
##Are social and caste toolkits enriched for developmental genes, as defined in Dmel?
##############
getOGG <- function(list, column){
  gene = rownames(list)[list$FDR < FDR]
  oggC = ogg111$OGG[ogg111[,column] %in% gene]
  return(oggC)
}

checkOver <- function(oggC,numC){
  MasterOGG = oggDevF$OGG[oggDevF$OGG %in% oggC]
  a = rbind(c(length(MasterOGG),nrow(oggDevF) - length(MasterOGG)),
            c(length(oggC) - length(MasterOGG),
              nrow(ogg111) - nrow(oggDevF) - length(oggC) + length(MasterOGG)))
  r = fisher.test(a,alternative="greater")
  return(c(length(oggC),length(MasterOGG),r$p.value))
}

res = list(list(list()))
tests = c("caste","social")
bigL = list(casteList,socList)
names(bigL) = tests
cols = c(2,3)
names(cols) = c("ant","bee")

for (test in tests){
  for (spec in c("bee","ant")){
    for (tissue in tissues){
      oggC = getOGG(bigL[[test]][[spec]][[tissue]],cols[spec])
      res[[test]][[spec]][[tissue]]=checkOver(oggC)
    }  
  }
}


#############
##Principle component analysis
#############
allF <- rbind(factorB,factorA)
allF$species = "bee"
allF$species[grepl("Ant",allF$sample)]="ant"
allF$stage_tissue = do.call(paste, c(allF[,c("stage","tissue")],list(sep='-')))
allF$bigStage = 'adult'
allF$bigStage[allF$stage==7]='pupa'
allF$bigStage[allF$stage!=7&allF$stage!=8]='brood'
allT = merge(beeTO,antTO,by="OGG")
rownames(allT) = allT$OGG
allT = allT[,-c(1)]

all_trans = log(allT + sqrt(allT ^ 2 + 1))
pc <- prcomp(t(all_trans),scale=T)
cor(as.numeric(as.factor(allF$species)),pc$x[,1])
cor(as.numeric(as.factor(allF$bigStage)),pc$x[,2])
cor(as.numeric(as.factor(allF$tissue)),pc$x[,7])

data = as.data.frame(pc$x)
data = cbind(data,allF)
data$stage_tissue = do.call(paste, c(allF[,c("stage","tissue")],list(sep='-')))
ggplot(data,aes(x=PC1,y=PC2,color=bigStage,shape=species))+
  geom_point()
  
ggplot(data,aes(x=PC2,y=PC3,color=stage,shape=tissue))+
  geom_point()

antT_l = log(antT+sqrt(antT^2+1))
pc <- prcomp(t(antT_l[rowSums(antT_l) > 0,]),scale=T)
data = as.data.frame(pc$x)
data = cbind(data,factorA)
ggplot(data[data$tissue=="larva",],aes(x=PC1,y=PC3,color=stage,shape=colony))+
  geom_point()
#1st component is species, then stage, then tissue




