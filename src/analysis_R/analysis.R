######
##Libraries and constants
######

library(edgeR)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(plyr)
library(ppcor)
library(viridis)
library(ggmosaic)

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

#Constants
FDR = 0.05
tissues <- c("head","gaster","mesosoma")

######
##Functions
######
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
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male"&!grepl("_V",f$sample),]) #remove eggs and 1st instars as well as males
  f$tissue_stage = as.factor(do.call(paste,c(f[,c(2,3)],list(sep="_"))))#Want to treat each larval stage separately, as well as each tissue of the adult stage. It's necessary to 
  #do this because the model can't estimate the effect of stage and tissue independently, since there are only adult h/m/g and only larvae of stage 3-6
  relev = levels(f$tissue_stage)
  relev[6:8]=c("8_head","8_mesosoma","8_gaster") #intially alphabetical
  f$tissue_stage = factor(f$tissue_stage,levels=relev)
  return(f)
}

#Treat larva as one stage
editFactor_oneLarv <- function(f){
  f <- droplevels(f[f$stage!=1&f$stage!=2&f$caste!="male"&!grepl("_V",f$sample),]) #remove eggs and 1st instars as well as males
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

#Return list of genes that are DE for caste overall
speciesCasteDE <- function(d,f,factorFun){
  f = factorFun(f)
  d <- d[,colnames(d) %in% f$sample]
  design <- model.matrix(~caste+tissue_stage+colony+caste*tissue_stage,data=f) #Note that we are lumping virgin queens, as well as nurses & foragers
  int <- EdgeR(d,design,(4+length(levels(f$tissue_stage))):(2+length(levels(f$tissue_stage))*2)) 
  genesInt <- rownames(int)[int$FDR < FDR]
  design <- model.matrix(~caste+tissue_stage+colony,data=f) 
  caste <- EdgeR(d,design,2) 
  genesDE <- rownames(caste)[caste$FDR < FDR]
  design <- model.matrix(~caste+tissue_stage+colony,data=f) 
  caste2 <- EdgeR(d[!rownames(d) %in% genesInt,],design,2)  #Remove genes with a significant interaction factor
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
    caste <- EdgeR(ds,design,2)
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

#Run all DE tests
speciesDEtests <- function(d,f,factorFun){
  OverallDE <- speciesCasteDE(d,f,factorFun)
  Stage <- speciesCasteStage(d,f,factorFun)
  ConsistentDE <- speciesCasteAllStage(Stage)
  return(list(OverallDE,Stage,ConsistentDE))
}

#Return logFC and FDR across development
collapseLogFC <- function(resList){
  d1 = d2 = data.frame(Gene = rownames(resList[[1]][[2]]))
  d1$L2 = resList[[1]][[2]]$logFC
  d2$L2 = resList[[1]][[2]]$FDR
  for (i in 2:length(resList)){
    dN = resList[[i]][[2]]
    dN$Gene = rownames(dN)
    colnames(dN)[c(1,5)] = rep(names(resList)[i],2)
    d1 = merge(d1,dN[,c(1,6)],by = "Gene")
    d2 = merge(d2,dN[,c(5,6)],by = "Gene")
  }
  return(list(d1,d2))
}

#Returns dataframe of DE direction based on fold change and FDR across development
DE_direction <- function(data){
  d = data.frame(Gene = data[[1]]$Gene)
  for (i in 2:ncol(data[[1]])){
    d[,i]  = "nonDE"
    d[,i][data[[1]][,i] > 0 & data[[2]][,i] < 0.05] = "worker"
    d[,i][data[[1]][,i] < 0 & data[[2]][,i] < 0.05] = "queen"
  }
  colnames(d) = colnames(data[[1]])
  return(d)
}

#Create heatmap of differential expression (number of times DE for queens and workers)
DEheatmap <- function(dfDE,species){
  dfDE$numQueen = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == "queen"))
  dfDE$numWorker = apply(dfDE[,c(2:ncol(dfDE))],1,function(x) sum(x == "worker"))
  d = table(dfDE$numQueen,dfDE$numWorker)
  m = melt(d)
  colnames(m)[c(1,2)] = c("QB","WB")
  p <- ggplot(m,aes(x=QB,y=WB))+
    geom_tile(aes(fill = value))+
    scale_fill_gradient(name = "count",trans = "log",
                        breaks = c(1,10,100,1000,10000),
                        labels = c(1,10,100,1000,10000))+
    geom_text(aes(x = QB,y = WB,label = value),color="white")+
    main_theme+
    scale_y_continuous(name = "number of times worker-biased",
                       breaks  = seq(0,max(m$WB)),
                       expand = c(0,0))+
    scale_x_continuous(name = "number of times queen-biased",
                       breaks = seq(0,max(m$QB)),
                       expand = c(0,0))+
    ggtitle(species)+
    theme(legend.position = "right",
          axis.line=element_line(color="black"),
          axis.text = element_text(size=16),
          axis.title = element_text(size = 22,face="bold"),
          plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
          panel.border = element_rect(size = 1, color = "black",fill = NA)) 
  return(p)
}

#Correlation of log fold change at each stage. First dataframe is signed, second is unsigned (i.e. absolute value)
lfcCor <- function(antD,beeD){
  nStage = ncol(antD) - 1 
  antD <- merge(antD,ogg11,by.x = "Gene",by.y = "gene_Mphar")
  beeD <- merge(beeD,ogg11,by.x = "Gene",by.y = "gene_Amel")
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

#Return list of nurse/forager DE genes for each tissue
speciesSocial <- function(d,f){
  f = droplevels(f[!is.na(f$NF),])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~NF+colony,data=fs)
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}

#Identify genes DE across development for the given subset of samples
casteDev <- function(fdr,caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(rownames(out)[out$FDR < fdr])
}

#Differential expression based on sex (male vs mated queen)
sexDE <- function(d,f){
  f = droplevels(f[grepl("AQ|_M",f$sample) & f$stage==8,])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~caste+colony,data=fs) 
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}

#Differential expression based on ovary activation (virgin vs mated queen)
ovaryDE <- function(d,f){
  f = droplevels(f[grepl("AQ|VQ",f$sample),])
  d <- d[,colnames(d) %in% f$sample]
  results = list()
  for (lev in levels(f$tissue)){
    fs = droplevels(f[grepl(as.character(lev),f$tissue),])
    ds <- d[,colnames(d) %in% rownames(fs)]
    design <- model.matrix(~VM+colony,data=fs) 
    results[[lev]] <- EdgeR(ds,design,2)
  }
  return(results)
}

#Generate developmental toolkit: genes DE across develoment in queens and workers
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

#Generate confidence intervals through bootstrapping
bootCI <- function(v,boots){
  v = v[!is.na(v)]
  bV = unlist(lapply(seq(1,boots),function(x){
    mean(v[sample(x = seq(1,length(v)),size = length(v),replace=TRUE)])
  }))
  c1 = quantile(bV, 0.025)
  c2 = quantile(bV, 0.975)
  m = mean(v)
  return(data.frame(mean = m, c1 = c1, c2 = c2))
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

#Input named (by gene) vector of coefficient of variation or another variable and output scatterplot
cb_cv <- function(data,cv,species){
  data$caste_bias = apply(data[,c(2:ncol(data))],1,function(x) sqrt(sum(x^2)))
  cv = as.data.frame(cv)
  cv$Gene = rownames(cv)
  data = merge(data,cv,by="Gene")
  p1 <- ggplot(data,aes(x = cv, y = caste_bias))+
    geom_point(alpha = 0.4)+
    geom_smooth(method = "loess",color = "red")+
    main_theme+
    xlab("mean coefficient of variation")+
    ylab("total caste bias")+
    ggtitle(species)
  return(list(p1,data))
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

#Get density of paired x and y values
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#Plot density of x and y values with trendline
density_plot <- function(stats,x,y){
  stats = stats[!is.na(stats[,x]) & !is.na(stats[,y]),]
  stats$density=get_density(stats[,x],stats[,y])
  p <- ggplot(stats,aes(x=stats[,x],y=stats[,y],color=density))+
    geom_point()+
    geom_smooth()+
    scale_color_viridis()+
    ylab(y)+
    xlab(x)+
    main_theme
  return(p)
}

#For Drosophila, Determine the correlation of log fold change and the distance between i and j
weightedLogFC <- function(gene){
  fc = c()
  dist = c()
  for (i in 1:length(gene)){
    for (j in 1:length(gene)){
      if (i >= j){
        next;
      }
      fc = c(fc,abs(log2(gene[i]/gene[j])))
      dist = c(dist,abs(i-j))
    } 
  }
  d = data.frame(fc = as.numeric(as.character(fc)), dist = as.numeric(as.character(dist)))
  d = d[!is.na(d$fc) & !is.infinite(d$fc),]
  if (nrow(d) < 2){
    return(NA)
  }
  test = cor.test(d$fc,d$dist,method="spearman")
  return(test$p.value)
}

#Caste-associated genes are also sex-biased
extractBias <- function(DEres){
  sexQ <- rownames(DEres)[DEres$FDR < 0.05 & DEres$logFC < 0]
  sexM <- rownames(DEres)[DEres$FDR < 0.05 & DEres$logFC > 0]
  sexFC <- data.frame(Gene = rownames(DEres), FC = DEres$logFC)
  return(list(FC = sexFC,Queen = sexQ,nonQueen = sexM))
}

#Generate plot of log fold change of two DE results, with DE genes highlighted
FCplot <- function(test1,test2){
  FC = merge(test1[[1]],test2[[1]],by = "Gene")
  FC$DE = "nonDE"
  FC$DE[FC$Gene %in% c(test1[[2]],test2[[2]],test1[[3]],test2[[3]])] = "DE-inconsistent"
  FC$DE[(FC$Gene %in% test1[[2]] & FC$Gene %in% test2[[3]]) | (FC$Gene %in% test1[[3]] & FC$Gene %in% test2[[2]])] = "DE-opp"
  FC$DE[FC$Gene %in% test1[[2]] &FC$Gene %in% test2[[2]]] = "Queen-upreg"
  FC$DE[FC$Gene %in% test1[[3]] &FC$Gene %in% test2[[3]]] = "Queen-downreg"
  p <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
    geom_point(aes(color = DE))+
    geom_smooth()+
    theme_bw()
  return(list(p,cor.test(FC$FC.x,FC$FC.y,method = "spearman"),table(FC$DE)))
}


############
##Loading in, cleaning data
###########
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
oggFilt = ogg2[ogg2$OGG %in% names(t),] #get 1-1 orthologs
b = table(oggFilt$gene_Amel)
b = b[b==1]
oggFilt = oggFilt[oggFilt$gene_Amel %in% names(b),]
a = table(oggFilt$gene_Mphar)
a = a[a==1]
ogg11 = oggFilt[oggFilt$gene_Mphar %in% names(a),]

#get orthologs with drosophila
ogg111 <- read.csv("~/GitHub/devnetwork/data/ThreeWayOGGMap.csv")
#Generated from Drosophila melanogaster gff file with the command grep "CDS.*FBpp.*gene=" GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff | sed 's/.*FLYBASE:\(FBpp.*\),/\1,/' | sed 's/,.*gene=/\t/' | sed 's/;.*//' | uniq> test
isoMap <- read.table("~/GitHub/devnetwork/data/Dmel_CDStoGene_key.txt")

ogg111 = merge(ogg111,isoMap,by.x="gene_Dmel_FLYBASE",by.y="V1")

t = table(ogg111$OGG)
t = t[t==1]
t2 = table(ogg111$gene_Mphar)
t2 = t2[t2==1]
t3 = table(ogg111$gene_Amel)
t3 = t3[t3==1]
t4 = table(ogg111$V2)
t4 = t4[t4==1]
ogg3f = ogg111[ogg111$OGG %in% names(t) & ogg111$gene_Mphar %in% names(t2) &
                 ogg111$gene_Amel %in% names(t3) & ogg111$V2 %in% names(t4),]

########
##Identify differentially expressed genes
########
beeTests = speciesDEtests(bee,factorB,editFactor)
beeTests_oneLarv = speciesDEtests(bee,factorB,editFactor_oneLarv)
antTests = speciesDEtests(ant,factorA,editFactor)
antTests_oneLarv = speciesDEtests(ant,factorA,editFactor_oneLarv)
ant_sexDE = sexDE(ant,factorA) #Comparing mated queens to males
bee_sexDE = sexDE(bee,factorB)
ant_VM = ovaryDE(ant,factorA) #Comparing mated to virgin queens
bee_VM = ovaryDE(bee,factorB)
beeSocial <- speciesSocial(bee,factorB) #comparing nurses to foragers
antSocial <- speciesSocial(ant,factorA)

#Extract logFC, FDR from DE results
antRes = collapseLogFC(antTests_oneLarv[[2]])
beeRes = collapseLogFC(beeTests_oneLarv[[2]])

#Make dataframes for differential expression at each stage
antDE_allstage = DE_direction(antRes)
beeDE_allstage = DE_direction(beeRes)
antDE = DE_direction(antRes)
beeDE = DE_direction(beeRes)

#Get genes with interaction effects
intAnt = antTests[[1]][[2]]
intBee = beeTests[[1]][[2]]
ogg11$intAnt = ogg11$intBee = 0
ogg11$intAnt[ogg11$gene_Mphar %in% intAnt] = 1
ogg11$intBee[ogg11$gene_Amel %in% intBee] = 1

#Extract abdomen genes for later use
ogg11$abdAnt = ogg11$abdBee = "nonDE"
ogg11$abdAnt[ogg11$gene_Mphar %in% antDE$Gene[antDE$`8_gaster`=="worker"]]="worker"
ogg11$abdAnt[ogg11$gene_Mphar %in% antDE$Gene[antDE$`8_gaster`=="queen"]]="queen"
ogg11$abdBee[ogg11$gene_Amel %in% beeDE$Gene[beeDE$`8_gaster`=="worker"]]="worker"
ogg11$abdBee[ogg11$gene_Amel %in% beeDE$Gene[beeDE$`8_gaster`=="queen"]]="queen"

#Figure 1a,b: Caste-bias changes across development
png("~/GitHub/devnetwork/figures/AntDEnums.png",width=2000,height=2000,res=300)
DEheatmap(antDE,"ant")
dev.off()

png("~/GitHub/devnetwork/figures/BeeDEnums.png",width=2000,height=2000,res=300)
DEheatmap(beeDE,"bee")
dev.off()

#Compare DE definition at each stage
aM = melt(antDE,id.vars = "Gene")
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE = sum(value!="nonDE"))
bM = melt(beeDE,id.vars = "Gene")
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE = sum(value!="nonDE"))
colnames(bM)[3] = "value_apis"
aM = merge(aM, ogg11,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM, ogg11,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGG","variable"))

allD = ddply(allM,~variable,summarize,
             NDE_Mphar = sum(value=="nonDE"),
             DE_Mphar = sum(value!="nonDE"),
             NDE_Apis = sum(value_apis=="nonDE"),
             DE_Apis = sum(value_apis!="nonDE"),
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

allD$LS_OGG_Mphar = allD$DE_Mphar - allD$DEboth
allD$LS_OGG_Apis = allD$DE_Apis - allD$DEboth
allD$LS_noOGG_Mphar = aD$DE - allD$DE_Mphar
allD$LS_noOGG_Apis = bD$DE - allD$DE_Apis
allD$NDE_Mphar = aD$NDE
allD$NDE_Apis = bD$NDE

allD = allD[,-c(3,5)]
allD[,1] = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen")
colnames(allD)[1] = "Stage"
allDm = melt(allD,id.vars = "Stage")
allDm$species = "Mphar"
allDm$species[grepl("Apis",allDm$variable)] ="Amel"
sub = allDm[allDm$variable=="DEboth",]
sub$species = "Amel"
allDm = rbind(allDm,sub)
allDm$variable = gsub("_Mphar","",allDm$variable)
allDm$variable = gsub("_Apis","",allDm$variable)
allDm = droplevels(allDm)
allDm$variable = factor(allDm$variable,rev(c("DEboth","LS_OGG","LS_noOGG","NDE")))
allDm$Stage = factor(allDm$Stage,levels = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen"))

png("~/GitHub/devnetwork/figures/AntDE_shared.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Mphar",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity")+
  main_theme+
  ylab("number of genes")+
  ggtitle("ant")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        plot.margin = margin(1.75,1.75,1.75,1.75,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/BeeDE_shared.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Amel",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity")+
  main_theme+
  ylab("number of genes")+
  ggtitle("bee")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.2,0.8),
        legend.title = element_blank(),
        plot.margin = margin(1.75,1.75,1.75,1.75,"cm"))
dev.off()

mat = matrix(ncol=3)
colnames(mat) = c("Stage","variable","species")
for (i in 1:nrow(allDm)){
  temp = data.frame(Stage = rep(allDm$Stage[i],allDm$value[i]),
                    variable = rep(allDm$variable[i],allDm$value[i]),
                    species = rep(allDm$species[i],allDm$value[i]))
  mat = rbind(mat,temp)
}
mat = mat[-c(1),]
mat$Stage = factor(mat$Stage,levels = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen"))
mat$variable = factor(mat$variable, levels = c("DEboth","LS_OGG","LS_noOGG","NDE"))

#Figure 1c. Number of shared DEs at each stage
png("~/GitHub/devnetwork/figures/numSharedDE_Apis.png",height=2000,width=2000,res=300)
ggplot(data=mat[mat$variable!="NDE",])+
  geom_mosaic(aes(x = product(variable,Stage),
                  fill = factor(variable)))+
  facet_grid(species~.)+
  ylab("proportion of DE genes")+
  xlab("Stage")+
  main_theme+
  theme(axis.text.x = element_text(angle = -25,hjust=0.1))+
  scale_x_productlist(labels = levels(mat$Stage),breaks = c(0,0.08,0.2,0.4,0.7))
dev.off()

#Association between interaction effect genes
x <- vennCounts(ogg11[,c("intBee","intAnt")])

#Figure 1d
png("~/GitHub/devnetwork/figures/InteractionOverlap.png")
vennDiagram(x,names = c("Apis","Mphar"))
dev.off()

chisq.test(rbind(c(868,1493),c(1026,3338)))

#Compare logFC across development between bees and ants
d = lfcCor(antRes[[1]],beeRes[[1]])
d[[1]]$Stage=d[[2]]$Stage = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen")
d[[1]]$Stage = d[[2]]$Stage = factor(d[[1]]$Stage,levels=d[[1]]$Stage)
d[[1]]$type = "signed"
d[[2]]$type = "unsigned"
df = rbind(d[[1]],d[[2]])

#Figure 1e- log fold change comparison between species
png("~/GitHub/devnetwork/figures/lfc_both.png",width=2000,height=2000,res=300)
ggplot(df,aes(x = Stage, y = cor,fill=type))+
  geom_bar(stat = "identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),width = 0.2,position=position_dodge(width=.9))+
  ylab("pearson correlation")+
  xlab("stage/tissue")+
  scale_fill_manual(values = c("grey20","grey60"),name = "correlation type")+
  main_theme+
  theme(legend.position = c(0.3,0.85),
        axis.text.x = element_text(angle=-45,hjust=0),
        plot.margin = margin(1.5,1.5,1.5,1.5,"cm"),
        legend.background = element_blank())
dev.off()

#Look at comparison of logFC in abdomens between species
antAbd = antRes[[1]][,c(1,6)]
beeAbd = beeRes[[1]][,c(1,6)]
antAbd$antDE = beeAbd$beeDE = "nonDE"
antAbd$antDE[antDE$`8_gaster`!="nonDE"]="DE"
beeAbd$beeDE[beeDE$`8_gaster`!="nonDE"]="DE"

antAbd = merge(antAbd,ogg11,by.x = "Gene",by.y = "gene_Mphar")
beeAbd = merge(beeAbd,ogg11,by.x = "Gene",by.y = "gene_Amel")

colnames(beeAbd)[2] = "abdomen_bee"
Abd = merge(antAbd,beeAbd,by = "OGG")
Abd = merge(Abd, ogg3f,by.x = "Gene.x",by.y="gene_Mphar",all.x = TRUE)
Abd$DE = "nonDE"
Abd$DE[Abd$antDE=="DE" & Abd$beeDE == "DE"]="bothDE"
Abd$DE[Abd$antDE=="DE" & Abd$beeDE != "DE"]="antDE"
Abd$DE[Abd$antDE!="DE" & Abd$beeDE == "DE"]="beeDE"

colnames(Abd)[c(3,7)] = c("ant_logFC","bee_logFC")

#Figure 1f. Log fold-change between ants and bees in queen abdomens
png("~/GitHub/devnetwork/figures/lfc_abd.png",width = 2000, height = 2000, res = 300) 
ggplot(Abd,aes(x = -ant_logFC, y = -bee_logFC))+ #We inverse the axes because we want queen genes upregulated
  geom_point(aes(color=DE),alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE,color = "black")+
  ylab("bee")+
  xlab("ant")+
  ggtitle("log2 fold-change\n(queen abdomen/worker abdomen)")+
  main_theme+
  theme(legend.position = c(0.15,0.85),
        legend.title = element_blank())
dev.off()

#Making table of commonly DE genes
queen = Abd[Abd$ant_logFC < 0 & Abd$bee_logFC < 0 & Abd$DE == "bothDE",]
worker = Abd[Abd$ant_logFC > 0 & Abd$bee_logFC > 0 & Abd$DE == "bothDE",]
queen_w = Abd[Abd$ant_logFC > 0 & Abd$bee_logFC < 0 & Abd$DE == "bothDE",]
worker_q = Abd[Abd$ant_logFC < 0 & Abd$bee_logFC > 0 & Abd$DE == "bothDE",]
q_ext = merge(queen,ext,by.x = "Gene.x",by.y="Gene")
w_ext = merge(worker,ext,by.x = "Gene.x",by.y="Gene")


########Part 2

#Development comparison between apis and mphar
fdr = 0.05
BeeDev <- genDevTool(fdr,factorB,bee)
AntDev <- genDevTool(fdr,factorA,ant)

ogg11$DevBee = ogg11$DevAnt = 0
ogg11$DevBee[ogg11$gene_Amel %in% BeeDev] = 1
ogg11$DevAnt[ogg11$gene_Mphar %in% AntDev] = 1

x <- vennCounts(ogg11[,c("DevBee","DevAnt")])

#Figure 2a
png("~/GitHub/devnetwork/figures/DevelomentalOverlap.png")
vennDiagram(x,names = c("Apis","Mphar"))
dev.off()

t = rbind(c(x[4,3],x[3,3]),c(x[2,3],x[1,3]))
chisq.test(t) #Highly significant

#Are developmental genes caste-associated?
ogg11$BeeQueen = ogg11$BeeWorker = ogg11$AntQueen = ogg11$AntWorker = ogg11$AntCaste = ogg11$BeeCaste = 0
ogg11$BeeQueen[ogg11$abdBee=="queen"]=1
ogg11$BeeWorker[ogg11$abdBee=="worker"]=1
ogg11$AntQueen[ogg11$abdAnt=="queen"]=1
ogg11$AntWorker[ogg11$abdAnt=="worker"]=1
ogg11$AntCaste[ogg11$abdAnt!="nonDE"]=1
ogg11$BeeCaste[ogg11$abdBee!="nonDE"]=1
ogg11$conInt = ogg11$conQueen = ogg11$conWorker = ogg11$conCaste = ogg11$conDev = 0
ogg11$conInt[ogg11$intBee + ogg11$intAnt == 2] = 1
ogg11$conQueen[ogg11$AntQueen + ogg11$BeeQueen == 2] = 1
ogg11$conWorker[ogg11$AntWorker + ogg11$BeeWorker == 2] = 1
ogg11$conCaste[ogg11$AntCaste + ogg11$BeeCaste == 2] = 1
ogg11$conDev[ogg11$DevAnt + ogg11$DevBee == 2] = 1

x <- vennCounts(ogg11[,c("conDev","conCaste")])

#Figure 2b. 
png("~/GitHub/devnetwork/figures/DevCaste.png")
vennDiagram(x,names = c("development","caste"))
dev.off()

x <- vennCounts(ogg11[,c("conDev","conInt")])

#Figure 2c. 
png("~/GitHub/devnetwork/figures/DevInt.png")
vennDiagram(x,names = c("development","casteXstage"))
dev.off()
t = rbind(c(x[4,3],x[3,3]),c(x[2,3],x[1,3]))
chisq.test(t) #Highly significant

AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[[2]][c(3:5)],function(x) extractBias(x[[2]]))
AvmRes <- lapply(ant_VM,extractBias)
AsocRes <- lapply(antSocial,extractBias)
AFC <- lapply(c(1:3),function(x) FCplot(AsexRes[[x]],AcasteRes[[x]]))
AFC_vm <- lapply(c(1:3),function(x) FCplot(AvmRes[[x]],AcasteRes[[x]]))
AFC_NF_caste <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AcasteRes[[x]]))
AFC_NF_sex <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AsexRes[[x]]))
AFC_NF_vm <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AvmRes[[x]]))

BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[[2]][c(3:5)],function(x) extractBias(x[[2]]))
BvmRes <- lapply(bee_VM,extractBias)
BsocRes <- lapply(beeSocial,extractBias)
BFC <- lapply(c(1:3),function(x) FCplot(BsexRes[[x]],BcasteRes[[x]]))
BFC_vm <- lapply(c(1:3),function(x) FCplot(BvmRes[[x]],BcasteRes[[x]]))
BFC_NF_caste <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BcasteRes[[x]]))
BFC_NF_sex <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BsexRes[[x]]))
BFC_NF_vm <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BvmRes[[x]]))


#Identify developmental genes in Drosophila based on RNA-seq across development 
counts <- read.csv("~/GitHub/devnetwork/data/counts_drosophila.csv")
fpkm <- read.csv("~/GitHub/devnetwork/data/fpkm_drosophila.csv")

factors <- read.csv("~/GitHub/devnetwork/data/drosophila_sra.csv")
rownames(counts) = rownames(fpkm) = counts$gene_id
counts = counts[,-c(1)]
fpkm = fpkm[,-c(1)]

keep = apply(fpkm,1,function(x) sum(x > 1) >= length(x)/2)
counts = counts[keep,]

factors$time_factor = as.factor(factors$time)
f_egg = factors[factors$time <= 12,]
f_devel = factors[factors$time <= 18,]
f_sex = factors[factors$time >= 25,]

eggCV = apply(counts[,colnames(counts) %in% f_egg$SRA],1,function(x) sd(x)/mean(x))
develCV = apply(counts[,colnames(counts) %in% f_devel$SRA],1,function(x) sd(x)/mean(x))

c_egg = counts[,f_egg$SRA]
c_egg = c_egg[,order(f_egg$time,decreasing=TRUE)] #Make samples ordered by time
c_devel = counts[,f_devel$SRA]
c_devel = c_devel[,order(f_devel$time,decreasing=TRUE)] #Make samples ordered by time

eggWLFC <- apply(c_egg,1,weightedLogFC) #P-value representing effect of development
develWLFC <- apply(c_devel,1,weightedLogFC)

#Sex-bias
design <- model.matrix(~sex+time_factor,data=droplevels(f_sex[f_sex$time>25,]))
d = counts[,colnames(counts) %in% f_sex$SRA[f_sex$time>25]]
sexGenes <- EdgeR(d,design,2)
sexGenes$Gene = rownames(sexGenes)

DmelDat = data.frame(Gene = names(eggWLFC),eggCV = eggCV, develCV = develCV, eggFC = eggWLFC, develFC = develWLFC)
DmelDat = merge(DmelDat,sexGenes,by = "Gene")

Dmel_oggDat = merge(DmelDat,ogg3f,by.x="Gene",by.y="V2")
Dogg = merge(Dmel_oggDat,ogg11,by="gene_Mphar")
DoggM = Dogg[,c(3:19,22:34)]
DoggM = melt(DoggM,id.vars = colnames(DoggM)[c(1:15)])
ggplot(DoggM,aes(x = variable, y = eggCV, fill = as.factor(value)))+
  geom_boxplot(notch=TRUE)

ggplot(DoggM,aes(x = variable, y = develCV, fill = as.factor(value)))+
  geom_boxplot(notch=TRUE)

ggplot(DoggM,aes(x = variable, y = -log(develFC), fill = as.factor(value)))+
  geom_boxplot(notch=TRUE)+coord_cartesian(ylim= c(0,5))


ggplot(DoggM,aes(x = variable, y = -log(eggFC), fill = as.factor(value)))+
  geom_boxplot(notch=TRUE)+coord_cartesian(ylim= c(0,5))

ggplot(DoggM,aes(x = variable, y = logFC, fill = as.factor(value)))+
  geom_boxplot(notch=TRUE)




bee_sub = bee[,grepl("_AQ|_N|_F",colnames(bee))]
pc <- prcomp(t(bee_sub))
group = rep("queen",ncol(bee_sub))
group[grepl("_N",colnames(bee_sub))]="nurse"
group[grepl("_F",colnames(bee_sub))]="forager"
tissue = rep("head",ncol(bee_sub))
tissue[grepl("M.",colnames(bee_sub))] = "mesosoma"
tissue[grepl("G.",colnames(bee_sub))] = "gaster"
d = data.frame(PC1 = pc$x[,1],PC2 = pc$x[,2],PC3 = pc$x[,3],group=group,tissue=tissue)
png("~/GitHub/devnetwork/figures/PCA_bee_adult.png",height=2000,width=2000,res=300)
ggplot(d,aes(x = PC1,y = PC2,color=group,shape=tissue))+
  geom_point()+theme_bw()
dev.off()


dmel_develGenes = names(eggWLFC)[eggWLFC < 0.2/length(eggWLFC)]

#Compare to hymenopteran developmental genes
devA = merge(oggDev,ogg3f,by = "gene_Mphar")
devA$DevDmel = 0
devA$DevDmel[devA$V2 %in% dmel_develGenes]=1

x <- vennCounts(devA[,c("DevApis","DevMphar","DevDmel")])

#Figure 2a
png("~/GitHub/devnetwork/figures/DevelomentalOverlap3.png")
vennDiagram(x,names = c("Apis","Mphar"))
dev.off()

fisher.test(rbind(c(51,483),c(127+81+39,nrow(devA) - 127 - 39 - 81 - 483 + 51)),alternative = "greater")

#Finding: There are no "consistent" DE
nGenes_consistentlyDE <- lapply(list(beeTests,beeTests_oneLarv,antTests,antTests_oneLarv),function(x){
  lapply(c(1,2,3),function(j) length(x[[3]][[j]]))
})

#There are, however some genes that are DE overall according to the model. We'll incorporate those below

############
###Part 3: Caste and social toolkits--in each species and the overlap
############


#########
##Calculating coefficient of variation
#########
factorA$meta = as.factor(apply(factorA[,c(2:6)],1,paste,collapse='_'))
factorB$meta = as.factor(apply(factorB[,c(2:6)],1,paste,collapse='_'))

antCV = averageCV(factorA,antT)
beeCV = averageCV(factorB,beeT)

p1 <- cb_cv(antRes[[1]],antCV[[1]],"ant")
p2 <- cb_cv(beeRes[[1]],beeCV[[1]],"bee")

png("~/GitHub/devnetwork/figures/caste_cv_ant.png",width=2000,height=2000,res=300)
p1[[1]]
dev.off()

png("~/GitHub/devnetwork/figures/caste_cv_bee.png",width=2000,height=2000,res=300)
p2[[1]]
dev.off()

p1 <- cb_cv(antSocRes[[1]],antCV[[1]],"ant")
p2 <- cb_cv(beeSocRes[[1]],beeCV[[1]],"bee")

png("~/GitHub/devnetwork/figures/social_cv_ant.png",width=2000,height=2000,res=300)
p1[[1]]+ylab("total nurse/forager bias")
dev.off()

png("~/GitHub/devnetwork/figures/social_cv_bee.png",width=2000,height=2000,res=300)
p2[[1]]+ylab("total nurse/forager bias")
dev.off()

p1 <- cb_cv(antRes[[1]][,c(1,7:9)],antCV[[1]],"ant")
p2 <- cb_cv(beeRes[[1]][,c(1,7:9)],beeCV[[1]],"bee")

cor.test(p1[[2]]$caste_bias,p1[[2]]$cv,method = "spearman")
cor.test(p2[[2]]$caste_bias,p2[[2]]$cv,method = "spearman")

#adding Jasper data
counts <- read.csv("~/GitHub/devnetwork/data/Jasper_counts.csv")
fpkm <- read.csv("~/GitHub/devnetwork/data/Jasper_fpkm.csv")
samples <- read.table("~/GitHub/devnetwork/data/samples_jasper.txt",head=T)
rownames(counts) = counts$gene_id
counts = counts[,-c(1)]
rownames(fpkm) = fpkm$gene_id
fpkm = fpkm[,-c(1)]

samples = samples[samples$caste!="callow",]
counts = counts[,colnames(counts) %in% samples$SRA]

tissueDE <- function(tissue){
  f = droplevels(samples[samples$tissue==tissue,])
  d = counts[,colnames(counts) %in% f$SRA]
  design <- model.matrix(~caste+rep,data = f)
  res = EdgeR(d,design,2)
  return(res)
}

brain <- tissueDE("brain")
antenna <- tissueDE("antenna")
hypophar <- tissueDE("hypophar")
mandibular <- tissueDE("mandibular")

mrjp <- lapply(list(brain,antenna,hypophar,mandibular),function(x) x[grepl("Mrjp",rownames(x)),])

samples$meta = as.factor(apply(samples[,c(1,2)],1,paste,collapse="_"))
colnames(samples)[4] = "sample"

bee_tissueCV = averageCV(samples,fpkm)



beeTissue <- calculatetau(samples,fpkm)
tau <- beeTissue[[1]]

p1 <- cb_cv(beeRes[[1]],tau,"bee")
png("~/GitHub/devnetwork/figures/tissueSpec_caste.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("tissue specificity (tau)")
dev.off()

p1 <- cb_cv(beeSocRes[[1]],tau,"bee")
png("~/GitHub/devnetwork/figures/tissueSpec_social.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("tissue specificity (tau)")+ylab("total behavior (nurse/forager) bias")
dev.off()

cor.test(p1[[2]]$caste_bias,p1[[2]]$cv,method = "spearman")

#Identify orthologs in bees to check for tissue-specificity 
AntOGGRes <- merge(antRes[[1]],ogg11,by.x="Gene",by.y="gene_Mphar")
AntOGGRes <- AntOGGRes[,-c(1,10)]
AntOGGRes <- AntOGGRes[,c(9,1:8)]
colnames(AntOGGRes)[1] = "Gene"
p1 <- cb_cv(AntOGGRes,tau,"ant")
png("~/GitHub/devnetwork/figures/tissueSpec_caste_ant.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("tissue specificity (tau)")+ylab("total caste bias")
dev.off()





devFCBee <- develIndex(bee,factorB)
devFCAnt <- develIndex(ant,factorA)

p1 <- cb_cv(antRes[[1]],devFCAnt,"ant")
png("~/GitHub/devnetwork/figures/develBias_caste_ant.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("total developmental bias")
dev.off()

cor.test(p1[[2]]$caste_bias,p1[[2]]$cv)

p1 <- cb_cv(beeRes[[1]],devFCBee,"bee")
png("~/GitHub/devnetwork/figures/develBias_caste_bee.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("total developmental bias")
dev.off()

p1 <- cb_cv(antSocRes[[1]],devFCAnt,"ant")
png("~/GitHub/devnetwork/figures/develBias_social_ant.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("total developmental bias")+ylab("total behavioral bias")
dev.off()

p1 <- cb_cv(beeSocRes[[1]],devFCBee,"bee")
png("~/GitHub/devnetwork/figures/develBias_social_bee.png",height=2000,width=2000,res=300)
p1[[1]]+xlab("total developmental bias")+ylab("total behavioral bias")
dev.off()









tauDF <- data.frame(Gene = names(tau),Tau=tau)
dA = antRes[[1]]
colnames(dA)[2:9] = c("L2","L3","L4","L5","Pupa","Adult_Head","Adult_Mesosoma","Adult_Abdomen")
dev = data.frame(Gene = names(devFCAnt), DevBias = devFCAnt)
cv = data.frame(Gene = names(antCV[[1]]),CV = antCV[[1]])
dA = merge(dA,dev,by="Gene")
dA = merge(dA,cv,by="Gene")
dA <- merge(dA,ogg11,by.x="Gene",by.y="gene_Mphar",all.x=TRUE)
dA <- merge(dA,tauDF,by.x="gene_Amel",by.y="Gene",all.x=TRUE)
dA$cb = apply(dA[,c(3:10)],1,function(x) sum(sqrt(x^2))) 
dA$cb_larv = apply(dA[,c(3:7)],1,function(x) sum(sqrt(x^2))) 
dA$cb_adult = apply(dA[,c(8:10)],1,function(x) sum(sqrt(x^2))) 
antExpr = data.frame(Gene = rownames(antT), Expr = apply(antT,1,function(x) log(mean(x)+1)))
dA = merge(dA,antExpr,by="Gene")
aStats = dA[,c(11,12,14:18)]
aPlots <- lapply(colnames(aStats[,-c(5,6)]), function(x){
  lapply(colnames(aStats[,-c(5,6)]), function(y) density_plot(aStats,x,y))
})

p <- unlist(aPlots,recursive=FALSE)
p <- lapply(p,function(x) x+theme(legend.position="none",
                                  axis.text=element_blank(),
                                  axis.title=element_blank()))
png("~/GitHub/devnetwork/figures/antAllCor.png",width=2000,height=2000,res=300)
do.call(grid.arrange,c(p,ncol=length(aPlots)))
dev.off()

cor(aStats,method="spearman")


dB = beeRes[[1]]
colnames(dB)[2:9] = c("L2","L3","L4","L5","Pupa","Adult_Head","Adult_Mesosoma","Adult_Abdomen")
dev = data.frame(Gene = names(devFCBee), DevBias = devFCBee)
cv = data.frame(Gene = names(beeCV[[1]]),CV = beeCV[[1]])
dB = merge(dB,dev,by="Gene")
dB = merge(dB,cv,by="Gene")
dB <- merge(dB,tauDF,by="Gene",all.x=TRUE)
dB$cb = apply(dB[,c(2:9)],1,function(x) sum(sqrt(x^2))) 
dB$cb_larv = apply(dB[,c(2:6)],1,function(x) sum(sqrt(x^2))) 
dB$cb_adult = apply(dB[,c(7:9)],1,function(x) sum(sqrt(x^2))) 
beeExpr = data.frame(Gene = rownames(beeT), Expr = apply(beeT,1,mean))
dB = merge(dB,beeExpr,by="Gene")
bStats = dB[,c(10:16)]
aPlots <- lapply(colnames(bStats), function(x){
  lapply(colnames(bStats), function(y) density_plot(aStats,x,y))
})

tissueExpr <- beeTissue[[2]]
tissueExpr$Gene = rownames(tissueExpr)
topTissue = apply(tissueExpr[,c(1:12)],1,function(x) colnames(tissueExpr)[1:12][x==max(x)])
topTissue <- lapply(topTissue,function(x){
  if (length(x) > 1) NA
  else x
})
t = unlist(topTissue)
tissueExpr$topTissue = t
dBTissue <- merge(dB,tissueExpr[,c(13,14)],by="Gene")
dBTissue = dBTissue[!is.na(dBTissue$Tau),]
dBTissue_topTau <- dBTissue[dBTissue$Tau > quantile(dBTissue$Tau,0.9),]
dBTissue_topTau$gland_type = "conserved"
dBTissue_topTau$gland_type[dBTissue_topTau$topTissue=="hypophar" | dBTissue_topTau$topTissue== "mandibular" |
                             dBTissue_topTau$topTissue == "nasonov" | dBTissue_topTau$topTissue == "sting" |
                             dBTissue_topTau$topTissue == "antenna"] = "novel"
dBTissue_topTau$gland_type = as.factor(dBTissue_topTau$gland_type)
png("~/GitHub/devnetwork/figures/topTissue_cb.png",height = 2000, width = 2000, res = 300)
ggplot(dBTissue_topTau,aes(x = topTissue,y=cb,fill = gland_type))+
  geom_boxplot()+
  main_theme+
  ylab("caste bias")+
  xlab("tissue with highest expression")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
dev.off()

png("~/GitHub/devnetwork/figures/topTissue_cb_cat.png",height = 2000, width = 2000, res = 300)
ggplot(dBTissue_topTau,aes(x = gland_type,y=cb,fill = gland_type))+
  geom_boxplot(notch=TRUE)+
  main_theme+
  theme(legend.position = "none")+
  ylab("caste bias")+
  xlab("tissue type")
dev.off()

wilcox.test(dBTissue_topTau$cb[dBTissue_topTau$gland_type=="novel"],dBTissue_topTau$cb[dBTissue_topTau$gland_type=="conserved"])

dBTissue <- merge(dA,tissueExpr[,c(13,14)],by.x = "gene_Amel",by.y="Gene")
dBTissue = dBTissue[!is.na(dBTissue$Tau),]
dBTissue_topTau <- dBTissue[dBTissue$Tau > quantile(dBTissue$Tau,0.9),]
dBTissue_topTau$gland_type = "conserved"
dBTissue_topTau$gland_type[dBTissue_topTau$topTissue=="hypophar" | dBTissue_topTau$topTissue== "mandibular" |
                             dBTissue_topTau$topTissue == "nasonov" | dBTissue_topTau$topTissue == "sting" |
                             dBTissue_topTau$topTissue == "antenna"] = "novel"
dBTissue_topTau$gland_type = as.factor(dBTissue_topTau$gland_type)
ggplot(dBTissue_topTau,aes(x = topTissue,y=cb,fill = gland_type))+geom_boxplot()
ggplot(dBTissue_topTau,aes(x = gland_type,y=cb,fill = gland_type))+geom_boxplot(notch=TRUE)
wilcox.test(dBTissue_topTau$cb[dBTissue_topTau$gland_type=="novel"],dBTissue_topTau$cb[dBTissue_topTau$gland_type=="conserved"])


p <- unlist(aPlots,recursive=FALSE)
p <- lapply(p,function(x) x+theme(legend.position="none",
                                  axis.text=element_blank(),
                                  axis.title=element_blank()))
png("~/GitHub/devnetwork/figures/beeAllCor.png",width=2000,height=2000,res=300)
do.call(grid.arrange,c(p,ncol=length(bStats)))
dev.off()


data = aStats
data = data[!is.na(data$Tau),]
pcor.test(data$cb,data$Tau,data$Expr,method="spearman")



##Sex-biased genes
design <- model.matrix(~sex+time_factor,data=droplevels(f_sex[f_sex$time>25,]))
d = counts[,colnames(counts) %in% f_sex$SRA[f_sex$time>25]]
sexGenes <- EdgeR(d,design,2)
sexDE <- rownames(sexGenes)[sexGenes$FDR < 0.001]
sexP <- data.frame(Gene = rownames(sexGenes), logP = -log(sexGenes$FDR))

DmelCV = data.frame(Gene = rownames(c_egg),eggCV = eggCV, develCV = develCV, egg_FC = eggWLFC, devel_FC = develWLFC)
DmelCV = merge(DmelCV,sexP,by="Gene")
dA_ogg <- merge(dA,ogg3f,by.x="Gene",by.y="gene_Mphar")
dA_ogg <- merge(dA_ogg,DmelCV,by.x="V2",by.y="Gene")
dB_ogg <- merge(dB,ogg3f,by.x="Gene",by.y="gene_Amel")
dB_ogg <- merge(dB_ogg,DmelCV,by.x="V2",by.y="Gene")

design <- model.matrix(~caste+colony,data = droplevels(factorA[grepl("AQH|_MH",factorA$sample),]))
d = ant[,colnames(ant) %in% rownames(design)]
antSex <- EdgeR(d,design,2)
aS = data.frame(Gene = rownames(antSex),logP_hymSex = -log(antSex$FDR))
dA_ogg = merge(dA_ogg,aS,by="Gene")

design <- model.matrix(~caste+colony,data = droplevels(factorB[grepl("AQH|_MH",factorA$sample),]))
d = bee[,colnames(bee) %in% rownames(design)]
beeSex <- EdgeR(d,design,2)
bS = data.frame(Gene = rownames(beeSex),logP_hymSex = -log(beeSex$FDR))
dB_ogg = merge(dB_ogg,bS,by="Gene")

corTable <- function(dA){
  xCols <- c("logP_hymSex","logP","eggCV","develCV","egg_FC","devel_FC")
  yCols <- c("DevBias","cb","cb_adult","cb_larv","logP_hymSex")
  results = matrix(nrow=6,ncol=5)
  results2 = matrix(nrow=6,ncol=5)
  
  for (i in 1:6){
    for (j in 1:5){
      if (i < 5){
        results[i,j] = cor.test(dA[,xCols[i]],dA[,yCols[j]],method="spearman",alternative="greater")$p.value
        results2[i,j] = cor.test(dA[,xCols[i]],dA[,yCols[j]],method="spearman",alternative="greater")$estimate
        
      } else {
        results[i,j] = cor.test(-log(dA[,xCols[i]]),dA[,yCols[j]],method="spearman",alternative="greater")$p.value
        results2[i,j] = cor.test(-log(dA[,xCols[i]]),dA[,yCols[j]],method="spearman",alternative="greater")$estimate
      }
    }
  }
  colnames(results) = colnames(results2) = yCols
  rownames(results) = rownames(results2) = xCols
  return(list(results,results2))
}

aTab = corTable(dA_ogg)
bTab = corTable(dB_ogg)

#Looking at development more discretely
ogg3f$ApisDev = ogg3f$MpharDev = ogg3f$DmelDev = ogg3f$hymDev = 0
ogg3f$ApisDev[ogg3f$gene_Amel %in% BeeDev] = 1
ogg3f$MpharDev[ogg3f$gene_Mphar %in% AntDev] = 1
ogg3f$DmelDev[ogg3f$V2 %in% DmelCV$Gene[DmelCV$egg_FC < 0.05/nrow(DmelCV)]]=1
ogg3f$hymDev[ogg3f$ApisDev+ogg3f$MpharDev==2]=1








