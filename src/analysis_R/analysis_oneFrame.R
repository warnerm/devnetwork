
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

#Identify genes DE across development for the given subset of samples
casteDev <- function(fdr,caste,data,factors){
  counts <- data[,(factors$tissue=="larva"|factors$tissue=="egg") & factors$caste==caste]
  f = droplevels(factors[factors$sample %in% colnames(counts),])
  design <- model.matrix(~stage+colony,data=f)
  out <- EdgeR(counts,design,2:length(levels(f$stage)))
  return(rownames(out)[out$FDR < fdr])
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
  FC$DE = "nonDE/inconsistent"
  #FC$DE[FC$Gene %in% c(test1[[2]],test2[[2]],test1[[3]],test2[[3]])] = "DE-inconsistent"
  #FC$DE[(FC$Gene %in% test1[[2]] & FC$Gene %in% test2[[3]]) | (FC$Gene %in% test1[[3]] & FC$Gene %in% test2[[2]])] = "DE-opp"
  FC$DE[FC$Gene %in% test1[[2]] &FC$Gene %in% test2[[2]]] = "Queen-upreg"
  FC$DE[FC$Gene %in% test1[[3]] &FC$Gene %in% test2[[3]]] = "Queen-downreg"
  p <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
    geom_point(aes(color = DE),alpha=0.5)+
    geom_smooth(method="lm",se=FALSE,color="black")+
    scale_color_manual(values = c("grey60","blue","red"))+
    main_theme
  return(list(p,cor.test(FC$FC.x,FC$FC.y,method = "spearman"),table(FC$DE)))
}


#Generate plot of log fold change of two DE results, with DE genes highlighted
FCplot_filt <- function(test1,test2){
  FC = merge(test1[[1]],test2[[1]],by = "Gene")
  FC$DE = "inconsistent"
  #FC$DE[FC$Gene %in% c(test1[[2]],test2[[2]],test1[[3]],test2[[3]])] = "DE-inconsistent"
  #FC$DE[(FC$Gene %in% test1[[2]] & FC$Gene %in% test2[[3]]) | (FC$Gene %in% test1[[3]] & FC$Gene %in% test2[[2]])] = "DE-opp"
  FC$DE[FC$Gene %in% test1[[2]] &FC$Gene %in% test2[[2]]] = "Queen-upreg"
  FC$DE[FC$Gene %in% test1[[3]] &FC$Gene %in% test2[[3]]] = "Queen-downreg"
  p <- ggplot(FC,aes(x = -FC.x,y = -FC.y))+ #Queen-up will be positive
    geom_point(aes(color = DE),alpha=0.5)+
    scale_color_manual(values = c("grey60","blue","red"))+
    main_theme
  return(list(p,cor.test(FC$FC.x,FC$FC.y,method = "spearman"),table(FC$DE)))
}


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


#Make venn diagram showing overlap of differential expression definitions
DEoverlap <- function(d, caste, species,larv = TRUE){
  if (larv){
    labels = c("L2","L3","L4","L5")
  } else {
    labels = c("larva","pupa","head","mesosoma","abdomen")
  }
  v = apply(d[,c(2:ncol(d))],2,function(x) d$Gene[x==caste])
  names(v) = labels
  venn.diagram(v,filename = paste("~/GitHub/devnetwork/figures/venn",caste,species,larv,".png",sep="_"),
               main = paste(species,caste,sep=", "),
               imagetype = "png",
               main.cex = 1.5)
}

DEoverlap_soc <- function(d, caste, species){
  labels = c("head","mesosoma","abdomen")
  v = apply(d[,c(2:ncol(d))],2,function(x) d$Gene[x==caste])
  names(v) = labels
  venn.diagram(v,filename = paste("~/GitHub/devnetwork/figures/venn",caste,species,".png",sep="_"),
               main = paste(species,caste,sep=", "),
               imagetype = "png",
               main.cex = 1.5)
}

mergeRes2 <- function(tax_id,species,taxNames,fullName,psFile){
  TGmap <- read.table(paste("~/GitHub/devnetwork/phylo_results/TGmap_",species,".txt",sep=""))
  blast <- read.table(paste("~/GitHub/devnetwork/phylo_results/",tax_id,"_",species,"_blastRes.tab",sep=""))
  PS <- read.table(paste("~/GitHub/devnetwork/phylo_results/",psFile,sep=""))
  ACUogg <- read.table(paste("~/GitHub/devnetwork/phylo_results/",species,"_acu",sep=""))
  
  mainPS = PS[!duplicated(PS),]
  oggBlast = merge(ACUogg,blast[,c(1,3)],by.x="V1",by.y="V3",all.x=TRUE)
  colnames(oggBlast)[3] = "transcript"
  oggPS = merge(oggBlast,mainPS,by.x="transcript",by.y="V1",all=TRUE)
  oggG = merge(oggPS,TGmap,by.x="transcript",by.y="V2")
  colnames(oggG) = c("transcript","ODBgene","OGGacu","ps","Gene")
  
  #Take oldest ps per gene
  oggGs <- ddply(oggG, ~ Gene + ODBgene + OGGacu,summarise,
                 ps = min(as.numeric(as.character(ps)),na.rm=TRUE))
  
  taxNames = c(taxNames,fullName)
  maxPS = length(taxNames)
  oggGs$ps[oggGs$ps==40]=maxPS
  oggGs$psName = apply(oggGs,1,function(x) taxNames[as.integer(as.character(x[4]))])
  
  return(oggGs)
}

#Plot fold change of phylostrata across development
pcFC <- function(res,ps,pLev){
  ps = ps[!is.na(ps$psName),]
  resM = melt(res,id.vars = "Gene")
  resM = resM[!is.na(resM$value),]
  resM = merge(ps,resM,by="Gene")
  resM$psName = factor(resM$psName,levels = pLev)
  p <- ggplot(resM,aes(x = variable,y=value,fill = psName))+
    geom_boxplot(outlier.shape=NA,notch = TRUE)+
    main_theme+
    geom_hline(yintercept = 0,color="black")+
    ylab("logFC (worker/queen)")+
    xlab("stage")
  return(p)
}

#Compare differential expression between species with fisher test
compDEfish <- function(antDE,beeDE){
  antDEo = merge(antDE,AB11,by.x="Gene",by.y="gene_Mphar")
  beeDEo = merge(beeDE,AB11,by.x="Gene",by.y="gene_Amel")
  results = list()
  for (col in colnames(antDEo[c(2:ncol(antDE))])){
    d1 = antDEo[,c(col,"OGGend")]
    d2 = beeDEo[,c(col,"OGGend")]
    d = merge(d1,d2,by='OGGend')
    colnames(d)[c(2,3)] = c("ant","bee")
    aDE = d$OGG[d$ant!="nonDE"]
    bDE = d$OGG[d$bee!="nonDE"]
    shared = sum(aDE %in% bDE)
    fishDE = fisher.test(rbind(c(shared,length(aDE) - shared),
                               c(length(bDE) - shared,nrow(d) - length(aDE) - length(bDE) + shared)),alternative = "greater")
    dShared = d[d$OGG %in% aDE[aDE %in% bDE],]
    aQ = dShared$OGG[dShared$ant=="queen"]
    bQ = dShared$OGG[dShared$bee=="queen"]
    bothQ = sum(aQ %in% bQ)
    fishCaste = fisher.test(rbind(c(bothQ,length(aQ) - bothQ),
                                  c(length(bQ) - bothQ,nrow(dShared) - length(aQ) - length(bQ) + bothQ)),alternative="greater")
    
    results[[col]] = c(shared,length(aDE),length(bDE),fishDE$estimate,fishDE$p.value,
                       bothQ,length(aQ),length(bQ),fishCaste$estimate,fishCaste$p.value)
  }
  
  res = ldply(results)
  colnames(res) = c("Stage/tissue","shared","antN","beeN",
                    "Fisher_odds","FisherP",
                    "sharedQ","antQ","beeQ","Fisher_odds","FisherP")
  
  res[,c(5,6,10,11)] = apply(res[,c(5,6,10,11)],2,function(x) signif(x,3))
  return(res)
}

#Get mean expression for each stage/caste combination
meanExpr <- function(d,f,factorFun){
  d = log(d^2+1)
  d = d[rowSums(d) > 0,]
  d = quantile_normalisation(d)
  f = editFactor(f)
  fQ = droplevels(f[f$caste=="queen",])
  fW = droplevels(f[f$caste=="worker",])
  fAll = list(fQ,fW)
  dAll = list(d[,colnames(d) %in% fQ$sample],d[,colnames(d) %in% fW$sample])
  expr = lapply(c(1,2),function(i){ lapply(levels(f$tissue_stage),function(x){
    rowSums(dAll[[i]][,colnames(dAll[[i]]) %in% fAll[[i]]$sample[fAll[[i]]$tissue_stage==x]])/sum(fAll[[i]]$tissue_stage==x)
  })})
  return(expr)
}

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

#Get pearson correlation of expression for a given sample
orthoExprCor <- function(exprA,exprB,N){
  dA = data.frame(Gene = names(exprA),exprAnt=exprA)
  dB = data.frame(Gene = names(exprB),exprBee = exprB)
  dA = merge(AB11,dA,by.x="gene_Mphar",by.y="Gene")
  dA = merge(dA,dB,by.x="gene_Amel",by.y="Gene")
  dA_order = dA[order(sqrt(dA$exprAnt*dA$exprBee),decreasing=TRUE),]
  dA_test = dA_order[1:N,]
  test = cor.test(dA_test$exprAnt,dA_test$exprBee)
  return(c(test$estimate,c1=test$conf.int[1],c2=test$conf.int[2]))
}

#Make plot of pearson correlation of expression for two lists of expression
pCorExpr <- function(mExprA,mExprB,N){
  results <- ldply(mapply(function(X,Y){
    ldply(lapply(seq(1,length(mExprA[[1]])),function(e){
        orthoExprCor(X[[e]],Y[[e]],N)
    }))
  },X=mExprA,Y=mExprB,SIMPLIFY = FALSE))
  results2 = ldply(lapply(seq(1,length(mExprA[[1]])),function(e){
    orthoExprCor(mExprA[[1]][[e]],mExprB[[2]][[e]],N)
  }))
  results3 = ldply(lapply(seq(1,length(mExprA[[1]])),function(e){
    orthoExprCor(mExprA[[2]][[e]],mExprB[[1]][[e]],N)
  }))
  results = do.call(rbind,list(results,results2,results3))
  results$caste = factor(rep(c("AQ-BQ","AW-BW","AQ-BW","AW-BQ"),each=length(mExprA[[1]])),levels=c("AQ-BQ","AW-BW","AQ-BW","AW-BQ"))
  results$stage = factor(rep(names(mExprA[[1]]),4),levels=names(mExprA[[1]]))
  p <- ggplot(results,aes(x = stage,y=cor,fill=caste))+
    geom_bar(stat="identity",position=position_dodge())+
    geom_errorbar(aes(ymin=c1,ymax=c2),width=0.3,position = position_dodge(width=0.9))+
    main_theme+
    ylab("pearson correlation")+
    ggtitle(paste("top",N,"genes",sep=" "))+
    theme(axis.text = element_text(angle=-25,hjust=0.1))
  return(p)
}

#merge endopterygota results
mergeEndo <- function(species,num){
  Dmel_END = read.table(paste("~/GitHub/devnetwork/phylo_results/",species,"_",num,"_END",sep=""))
  Dmel_blast =  read.table(paste("~/GitHub/devnetwork/phylo_results/",num,"_",species,"_blastRes.tab",sep=""))
  Dmel_tg = read.table(paste("~/GitHub/devnetwork/phylo_results/TGmap_",species,".txt",sep=""))
  
  Dmel_genes = merge(Dmel_blast,Dmel_tg,by.x="V1",by.y="V2")
  Dmel_genes = Dmel_genes[!duplicated(Dmel_genes[,4]),]
  Dmel_genes = merge(Dmel_genes,Dmel_END,by.x="V3",by.y="V1")
  DmelT = table(Dmel_genes$V2.y)
  Dmel_genes = Dmel_genes[Dmel_genes$V2.y %in% names(DmelT)[DmelT==1],][,c(4,5)]
  colnames(Dmel_genes) = c("Gene","OGGend")
  return(Dmel_genes)
  
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
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]

antT_filt = antT[rownames(antT) %in% rownames(ant),]
beeT_filt = beeT[rownames(beeT) %in% rownames(bee),]

write.csv(antT_filt,"~/GitHub/devnetwork/data/ant_filt.tpm.csv")
write.csv(beeT_filt,"~/GitHub/devnetwork/data/bee_filt.tpm.csv")


Atax = "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Holometabola; Hymenoptera; Apocrita; Aculeata; Formicoidea; Formicidae; Myrmicinae; Solenopsidini; Monomorium"
Atax = strsplit(Atax,"; ")[[1]]
Aps <- mergeRes2(307658,"Mphar",Atax,"novel","Mphar_307658_blastAll_ps")
ant_name = c("Solenopsidini","Formicidae","Myrmicinae")
old = c("Bilateria","Eumetazoa","Protostomia","Ecdysozoa","Metazoa")
insect_arthropod = c("Holometabola","Arthropoda","Hexapoda","Pancrustacea","Mandibulata","Pterygota","Neoptera")
aculeata = c("Aculeata","Apocrita")
bee_name = c("Apoidea", "Apidae", "Apinae","Apini")
Aps$psName[Aps$psName %in% ant_name] = "ant"
Aps$psName[Aps$psName %in% old] = "old"
Aps$psName[Aps$psName %in% insect_arthropod] = "insect_arthropod"
Aps$psName[Aps$psName %in% aculeata] = "aculeata"
Aps$psName[Aps$psName=="Monomorium"] = "novel"

Btax = "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Holometabola; Hymenoptera; Apocrita; Aculeata; Apoidea; Apidae; Apinae; Apini; Apis"
Btax = strsplit(Btax,"; ")[[1]]
Bps <- mergeRes2(7460,"Amel",Btax,"novel","Amel_7460_blastAll_ps")
Bps$psName[Bps$psName %in% old] = "old"
Bps$psName[Bps$psName %in% insect_arthropod] = "insect_arthropod"
Bps$psName[Bps$psName %in% aculeata] = "aculeata"
Bps$psName[Bps$psName %in% bee_name] = "bee"
Bps$psName[Bps$psName == "Apis"] = "novel"

AllPS = merge(Aps[!is.na(Aps$OGGacu),],Bps[!is.na(Bps$OGGacu),],by = "OGGacu")
AllPS1 = AllPS[,c(1,4,8)]
AllPS_sum = ddply(AllPS1,~OGGacu,summarize,
                  ps = min(ps.x,ps.y))
AllPS_sum$psName = "old"
AllPS_sum$psName[AllPS$ps > 8] = "insect"
AllPS_sum$psName[AllPS$ps > 18] = "hymenoptera"
AllPS_sum$psName[AllPS$ps > 20] = "aculeata"

AmelEND = mergeEndo("Amel",7460)
MpharEND = mergeEndo("Mphar",307658)
DmelEND = mergeEndo("Dmel",7227)
ENDogg = merge(AmelEND,MpharEND,by="OGGend")
ENDogg = merge(ENDogg,DmelEND,by="OGGend")
colnames(ENDogg) = c("OGGend","gene_Amel","gene_Mphar","gene_Dmel")

t = table(AllPS$OGGacu)
t = t[t==1]
AB11 = AllPS[AllPS$OGGacu %in% names(t),]
t = table(AB11$Gene.x)
t = t[t==1]
AB11 = AB11[AB11$Gene.x %in% names(t),]
t = table(AB11$Gene.y)
t = t[t==1]
AB11 = AB11[AB11$Gene.y %in% names(t),]
ACUogg = AB11[,c(1,2,6)]
colnames(ACUogg) = c("OGGacu","gene_Mphar","gene_Amel")

load("~/GitHub/devnetwork/data/DEtests.RData")

#Plot fold change of phylostrata across development

pcFC <- function(res,ps,pLev){
  ps = ps[!is.na(ps$psName),]
  resM = melt(res,id.vars = "Gene")
  resM = resM[!is.na(resM$value),]
  resM = merge(ps,resM,by="Gene")
  resM$psName = factor(resM$psName,levels = pLev)
  p <- ggplot(resM,aes(x = variable,y=value,fill = psName))+
    geom_boxplot(outlier.shape=NA,notch = TRUE)+
    main_theme+
    geom_hline(yintercept = 0,color="black")+
    ylab("logFC (worker/queen)")+
    xlab("stage")
  return(p)
}

pcFC_abd <- function(res,ps,pLev){
  ps = ps[!is.na(ps$psName),]
  resM = melt(res,id.vars = "Gene")
  resM = resM[!is.na(resM$value),]
  resM = merge(ps,resM,by="Gene")
  resM$psName = factor(resM$psName,levels = pLev)
  p <- ggplot(resM[resM$variable=="abdomen",],aes(x = psName,y=value))+
    geom_boxplot(outlier.shape=NA,notch = TRUE)+
    main_theme+
    geom_hline(yintercept = 0,color="black")+
    ylab("logFC (worker/queen)")+
    xlab("stage")
  return(p)
}

antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]])

antPS = merge(antRes[[1]],Aps,by="Gene")
antPS = merge(antPS,antCB,by="Gene")
beePS = merge(beeRes[[1]],Bps,by="Gene")
beePS = merge(beePS,beeCB,by="Gene")
antPS$abd_adj = antPS$abdomen - median(antPS$abdomen)
beePS$abd_adj = beePS$abdomen - median(beePS$abdomen)
antPS$species ="ant"
beePS$species = "bee"
allPS = merge(antPS,beePS,by="OGGacu")
allPS2 = rbind(antPS,beePS)
allPS2 = allPS2[!is.na(allPS2$psName),]
allPS2$psName = factor(allPS2$psName,levels = c("old","insect_arthropod","Hymenoptera","aculeata","ant",
                                                "bee","novel"))
p <- ggplot(allPS2,aes(x = species,y=abd_adj,fill = psName))+
  geom_boxplot(outlier.shape=NA,notch = TRUE)+
  main_theme+
  geom_hline(yintercept = 0,color="black")+
  ylab("logFC (worker/queen)")+
  xlab("stage")


png("~/GitHub/devnetwork/figures/psFC_bothabd.png",height=2000,width=2500,res=300)
p+coord_cartesian(ylim=c(-7,7))
dev.off()


antPS = merge(antRes[[2]],Aps,by="Gene")
antPS = merge(antPS,antCB,by="Gene")
beePS = merge(beeRes[[2]],Bps,by="Gene")
beePS = merge(beePS,beeCB,by="Gene")
aM = melt(antPS[,c(1:9,13,15,16,18)],id.vars = c("Gene","cb","cb_noAdult","psName","OGGacu"))
bM = melt(beePS[,c(1:9,13,15,16,18)],id.vars = c("Gene","cb","cb_noAdult","psName","OGGacu"))

ggplot(bM[bM$variable=="8_gaster",],aes(x = value,fill = psName))+
  geom_bar(stat="count",position = "fill")+
  main_theme+
  ylab("overall caste bias")+
  xlab("stage")


png("~/GitHub/devnetwork/figures/psFC_overall.png",height=2000,width=2000,res=300)
p
dev.off()


ApLev = c("old","insect_arthropod","Hymenoptera","aculeata","ant","novel")
png("~/GitHub/devnetwork/figures/psFC_ant.png",height=2000,width=2000,res=300)
pcFC(antRes[[1]],Aps,ApLev)+ggtitle("ant")+ylim(-6,6)+theme(legend.position = "top",legend.title = element_blank())
dev.off()

png("~/GitHub/devnetwork/figures/psFC_ant_abd.png",height=2000,width=2000,res=300)
pcFC_abd(antRes[[1]],Aps,ApLev)+ggtitle("ant")+ylim(-6,6)+theme(legend.position = "top",legend.title = element_blank())
dev.off()

png("~/GitHub/devnetwork/figures/psFC_ant_behav.png",height=2000,width=2000,res=300)
pcFC(antSocRes[[1]],Aps,ApLev)+ggtitle("ant behavior")+ylab("logFC (forager/nurse)")+theme(legend.position = "top",legend.title = element_blank())
dev.off()

BpLev = c("old","insect_arthropod","Hymenoptera","aculeata","bee","novel")
png("~/GitHub/devnetwork/figures/psFC_bee.png",height=2000,width=2000,res=300)
pcFC(beeRes[[1]],Bps,BpLev)+ggtitle("bee")+ylim(-6,6)+theme(legend.position = "top",legend.title=element_blank())
dev.off()

png("~/GitHub/devnetwork/figures/psFC_bee_abd.png",height=2000,width=2000,res=300)
pcFC_abd(beeRes[[1]],Bps,BpLev)+ggtitle("bee")+ylim(-6,6)+theme(legend.position = "top",legend.title=element_blank())
dev.off()


png("~/GitHub/devnetwork/figures/psFC_bee_behav.png",height=2000,width=2000,res=300)
pcFC(beeSocRes[[1]],Bps,BpLev)+ggtitle("bee behavior")+ylab("logFC (for/nurse)")+theme(legend.position = "top",legend.title=element_blank())
dev.off()

res = compDEfish(antRes[[2]],beeRes[[2]])
res[,1] = c("larva","pupa","adult_head","adult_thorax","adult_abdomen")

png("~/GitHub/devnetwork/figures/DEtable.png",height=2000,width=4000,res=300)
grid.table(res)
dev.off()

res = compDEfish(antSocRes[[2]],beeSocRes[[2]])
res[,1] = c("adult_head","adult_thorax","adult_abdomen")

png("~/GitHub/devnetwork/figures/DEtable_behav.png",height=2000,width=4000,res=300)
grid.table(res)
dev.off()

#Compare DE definition at each stage
antRes[[2]]$ortholog_found = antRes[[2]]$OGG_found = FALSE
antRes[[2]]$ortholog_found[antRes[[2]]$Gene %in% AllPS$Gene.x] = TRUE
antRes[[2]]$OGG_found[antRes[[2]]$Gene %in% ACUogg$gene_Mphar] = TRUE
aM = melt(antRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE_noOrtholog = sum(value!="nonDE" & !ortholog_found),
           DE_duplicate = sum(value!="nonDE" & ortholog_found & !OGG_found),
           DE_OGG = sum(value!="nonDE" & OGG_found))
beeRes[[2]]$ortholog_found = beeRes[[2]]$OGG_found = FALSE
beeRes[[2]]$ortholog_found[beeRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeRes[[2]]$OGG_found[beeRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE_noOrtholog = sum(value!="nonDE" & !ortholog_found),
           DE_duplicate = sum(value!="nonDE" & ortholog_found & !OGG_found),
           DE_OGG = sum(value!="nonDE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

aD$DEboth = bD$DEboth = allD$DEboth
aD$DE_OGG = aD$DE_OGG - aD$DEboth
bD$DE_OGG = bD$DE_OGG - allD$DEboth

aDM = melt(aD,id.vars = "variable")
bDM = melt(bD,id.vars = "variable")
colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")

png("~/GitHub/devnetwork/figures/AntDE_shared.png",height=2000,width=2000,res=300)
ggplot(aDM[aDM$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  coord_flip()+
  ylab("number of genes")+
  ggtitle("ant")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/BeeDE_shared.png",height=2000,width=2000,res=300)
ggplot(bDM[bDM$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  coord_flip()+
  ylab("number of genes")+
  ggtitle("bee")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/AntDE_shared_prop.png",height=2000,width=2000,res=300)
ggplot(aDM[aDM$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position="fill")+
  main_theme+
  coord_flip()+
  ylab("proportion of genes")+
  ggtitle("ant")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/BeeDE_shared_prop.png",height=2000,width=2000,res=300)
ggplot(bDM[bDM$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position="fill")+
  main_theme+
  coord_flip()+
  ylab("proportion of genes")+
  ggtitle("bee")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

allM_abd = allM[allM$variable == "abdomen",]
allM_ps = merge(allM_abd,AllPS,by="OGGacu")
allM_ps$DE = "nonDE"
allM_ps$DE[allM_ps$value!="nonDE" | allM_ps$value != "nonDE"] = "inconsistent"
allM_ps$DE[allM_ps$value == "worker" & allM_ps$value_apis == "worker"] = "worker"
allM_ps$DE[allM_ps$value == "queen" & allM_ps$value_apis == "queen"] = "queen"

ggplot(allM_ps,aes(x = psName,fill=DE))+
  geom_bar(stat="count")

##Same as above for social
aM = melt(antSocRes[[2]],id.vars = "Gene")
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE = sum(value!="nonDE"))
bM = melt(beeSocRes[[2]],id.vars = "Gene")
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           DE = sum(value!="nonDE"))
colnames(bM)[3] = "value_apis"
aM = merge(aM, AB11,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM, AB11,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGend","variable"))

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
allD[,1] = c("adult_head","adult_thorax","adult_abdomen")
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
allDm$Stage = factor(allDm$Stage,levels = c("larva","pupa","adult_head","adult_thorax","adult_abdomen"))

png("~/GitHub/devnetwork/figures/AntDE_shared_soc.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Mphar",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity")+
  main_theme+
  coord_flip()+
  ylab("number of genes")+
  ggtitle("ant")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/BeeDE_shared_soc.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Amel",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity")+
  main_theme+
  coord_flip()+
  ylab("number of genes")+
  ggtitle("bee")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/AntDE_shared_prop_soc.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Mphar",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity",position="fill")+
  main_theme+
  coord_flip()+
  ylab("proportion of genes")+
  ggtitle("ant")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

png("~/GitHub/devnetwork/figures/BeeDE_shared_prop_soc.png",height=2000,width=2000,res=300)
ggplot(allDm[allDm$variable!="NDE" & allDm$species=="Amel",],
       aes(x = Stage, y = value, fill = variable))+
  geom_bar(stat="identity",position="fill")+
  main_theme+
  coord_flip()+
  ylab("proportion of genes")+
  ggtitle("bee")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = c(0.7,0.2),
        legend.title = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"))
dev.off()

#Identify genes with conserved caste bias in abdomens
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

allMp <- merge(allM,Aps,by.x="gene_Mphar",by.y="Gene")

allD = ddply(allMp,~variable+psName,summarize,
             NDE_Mphar = sum(value=="nonDE"),
             DE_Mphar = sum(value!="nonDE"),
             NDE_Apis = sum(value_apis=="nonDE"),
             DE_Apis = sum(value_apis!="nonDE"),
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))


#Compare logFC across development between bees and ants
d = lfcCor(antRes_allstage[[1]],beeRes_allstage[[1]])
#d[[1]]$Stage=d[[2]]$Stage = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen")
d[[1]]$Stage=d[[2]]$Stage = c("L2","L3","L4","L5","pupa","adult_head","adult_mesosoma","adult_abdomen")
d[[1]]$Stage = d[[2]]$Stage = factor(d[[1]]$Stage,levels=d[[1]]$Stage)
d[[1]]$type = "signed"
d[[2]]$type = "unsigned"
#df = rbind(d[[1]],d[[2]])

#Figure 1e- log fold change comparison between species
png("~/GitHub/devnetwork/figures/lfc_both.png",width=2000,height=2000,res=300)
ggplot(d[[1]],aes(x = Stage, y = cor))+
  geom_bar(stat = "identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),width = 0.2,position=position_dodge(width=.9))+
  ylab("pearson correlation")+
  xlab("stage/tissue")+
  main_theme+
  theme(legend.position = c(0.3,0.85),
        axis.text.x = element_text(angle=-45,hjust=0),
        plot.margin = margin(1.5,2,1.5,1.5,"cm"),
        legend.background = element_blank())
dev.off()

antRes_allstage[[3]][,c(2:9)] = -log(antRes_allstage[[3]][,c(2:9)]+1)
beeRes_allstage[[3]][,c(2:9)] = -log(beeRes_allstage[[3]][,c(2:9)]+1)
d = lfcCor(antRes_allstage[[3]],beeRes_allstage[[3]])
#d[[1]]$Stage=d[[2]]$Stage = c("larva","pupa","adult_head","adult_mesosoma","adult_abdomen")
d[[1]]$Stage=d[[2]]$Stage = c("L2","L3","L4","L5","pupa","adult_head","adult_mesosoma","adult_abdomen")
d[[1]]$Stage = d[[2]]$Stage = factor(d[[1]]$Stage,levels=d[[1]]$Stage)
d[[1]]$type = "signed"
d[[2]]$type = "unsigned"
#df = rbind(d[[1]],d[[2]])

#Figure 1e- log fold change comparison between species
png("~/GitHub/devnetwork/figures/lfc_both.png",width=2000,height=2000,res=300)
ggplot(d[[1]],aes(x = Stage, y = cor))+
  geom_bar(stat = "identity",position=position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),width = 0.2,position=position_dodge(width=.9))+
  ylab("pearson correlation")+
  xlab("stage/tissue")+
  main_theme+
  theme(legend.position = c(0.3,0.85),
        axis.text.x = element_text(angle=-45,hjust=0),
        plot.margin = margin(1.5,2,1.5,1.5,"cm"),
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


allPS <- merge(Aps[!is.na(Aps$OGGend),],Bps[!is.na(Aps$OGGend),],by="OGGend")

######Comparing to sex-biased and virgin-mated biased genes

#Three-way venn for each species


#log FC sex vs caste, v-m vs caste

#Are conserved caste-biased genes more likely to be sex-biased? Could look at correlation just for conserved genes




########Part 2

#Development comparison between apis and mphar
fdr = 0.05
BeeDev <- genDevTool(fdr,factorB,bee)
AntDev <- genDevTool(fdr,factorA,ant)
antDE = antRes[[2]]
beeDE = beeRes[[2]]
beeDE$devel = antDE$devel = "nonDE"
beeDE$devel[beeDE$Gene %in% BeeDev] = "DE"
antDE$devel[antDE$Gene %in% AntDev] = "DE"
antDE$devel=factor(antDE$devel,levels = c("nonDE","DE"))
beeDE$devel=factor(beeDE$devel,levels = c("nonDE","DE"))
aM = melt(antDE,id.vars = c("Gene","devel"))
aM$species = "ant"
bM = melt(beeDE,id.vars = c("Gene","devel"))
bM$species = "bee"

overlapPlot <- function(data){
  p <- ggplot(data,aes(x = value,fill=devel))+
    geom_bar(stat="count",position="fill")+
    facet_wrap(~variable,nrow=1)+
    main_theme
  return(p)
}

png("~/GitHub/devnetwork/figures/DevDE_caste_ant.png",width=4000,height=2000,res=300)
overlapPlot(aM)+ylab("proportion")
dev.off()


png("~/GitHub/devnetwork/figures/DevDE_caste_bee.png",width=4000,height=2000,res=300)
overlapPlot(bM)+ylab("proportion")
dev.off()


name1 = "queen"
name2 = "worker"
Mall = merge(aM,AB11,by.x="Gene",by.y="gene_Mphar")
Mall2 = merge(bM,AB11,by.x="Gene",by.y="gene_Amel")
Mall = rbind(Mall[,-c(7)],Mall2[,-c(7)])
Ma = merge(Mall2,aM,by.x=c("gene_Mphar","variable"),by.y=c("Gene","variable"))
Ma$cons = "nonDE"
Ma$cons[Ma$value.x==name1 & Ma$value.y==name1] = name1
Ma$cons[Ma$value.x==name2 & Ma$value.y==name2] = name2
Ma$cons[(Ma$value.x==name2 & Ma$value.y==name1) | (Ma$value.x==name1 & Ma$value.y==name2)] = "inconsistent"
Ma$devel_cons = "nonDE"
Ma$devel_cons[Ma$devel.x!= "nonDE" | Ma$devel.y != "nonDE"] = "inconsistent"
Ma$devel_cons[Ma$devel.x!= "nonDE" & Ma$devel.y != "nonDE"] = "DE"
Ma2 = Ma[,c("variable","OGGend","cons","devel_cons")]
colnames(Ma2)[c(3,4)] = c("value","devel")
Ma2$value = factor(Ma2$value,levels = c("nonDE","inconsistent",name1,name2))
Ma2$devel = factor(Ma2$devel,levels = c("nonDE","inconsistent","DE"))

png("~/GitHub/devnetwork/figures/DevDE_caste_conserved.png",width=4000,height=2000,res=300)
overlapPlot(Ma2[Ma2$variable!="larva" & Ma2$variable!="pupa",])+ylab("proportion")
dev.off()


antDE = antSocRes[[2]]
beeDE = beeSocRes[[2]]
beeDE$devel = antDE$devel = "nonDE"
beeDE$devel[beeDE$Gene %in% BeeDev] = "DE"
antDE$devel[antDE$Gene %in% AntDev] = "DE"
antDE$devel=factor(antDE$devel,levels = c("nonDE","DE"))
beeDE$devel=factor(beeDE$devel,levels = c("nonDE","DE"))
aM = melt(antDE,id.vars = c("Gene","devel"))
aM$species = "ant"
bM = melt(beeDE,id.vars = c("Gene","devel"))
bM$species = "bee"
aM$value = factor(aM$value,levels=c("nonDE","nurse","forager"))
bM$value = factor(bM$value,levels=c("nonDE","nurse","forager"))

png("~/GitHub/devnetwork/figures/DevDE_social_ant.png",width=4000,height=2000,res=300)
overlapPlot(aM)+ylab("proportion")
dev.off()


png("~/GitHub/devnetwork/figures/DevDE_social_bee.png",width=4000,height=2000,res=300)
overlapPlot(bM)+ylab("proportion")
dev.off()

name1 = "nurse"
name2 = "forager"
Mall = merge(aM,AB11,by.x="Gene",by.y="gene_Mphar")
Mall2 = merge(bM,AB11,by.x="Gene",by.y="gene_Amel")
Mall = rbind(Mall[,-c(7)],Mall2[,-c(7)])
Ma = merge(Mall2,aM,by.x=c("gene_Mphar","variable"),by.y=c("Gene","variable"))
Ma$cons = "nonDE"
Ma$cons[Ma$value.x==name1 & Ma$value.y==name1] = name1
Ma$cons[Ma$value.x==name2 & Ma$value.y==name2] = name2
Ma$cons[(Ma$value.x==name2 & Ma$value.y==name1) | (Ma$value.x==name1 & Ma$value.y==name2)] = "inconsistent"
Ma$devel_cons = "nonDE"
Ma$devel_cons[Ma$devel.x!= "nonDE" | Ma$devel.y != "nonDE"] = "inconsistent"
Ma$devel_cons[Ma$devel.x!= "nonDE" & Ma$devel.y != "nonDE"] = "DE"
Ma2 = Ma[,c("variable","OGGend","cons","devel_cons")]
colnames(Ma2)[c(3,4)] = c("value","devel")
Ma2$value = factor(Ma2$value,levels = c("nonDE","inconsistent",name1,name2))
Ma2$devel = factor(Ma2$devel,levels = c("nonDE","inconsistent","DE"))

overlapPlot(Ma2[Ma2$variable!="larva" & Ma2$variable!="pupa",])

MallD = ddply(Mall, ~ OGGend + variable,summarize,
              consWorker = sum(value=="forager") == 2,
              consQueen = sum(value=="nurse") == 2,
              consDevel = sum(devel=="DE") == 2)

p1 <- ggplot(data=MallD[MallD$consQueen==TRUE,])+
  geom_mosaic(aes(x = product(consDevel,variable),
                  fill = factor(consDevel)))+
  ylab("proportion of DE genes")+
  xlab("Stage")+
  main_theme+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = -25,hjust=0.1))

ogg11 = AB11
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
#Extract abdomen genes for later use
ogg11 = AB11
ogg11$DevBee = ogg11$DevAnt = 0
ogg11$DevBee[ogg11$gene_Amel %in% BeeDev] = 1
ogg11$DevAnt[ogg11$gene_Mphar %in% AntDev] = 1
ogg11$abdAnt = ogg11$abdBee = "nonDE"
ogg11$abdAnt[ogg11$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen =="worker"]]="worker"
ogg11$abdAnt[ogg11$gene_Mphar %in% antRes[[2]]$Gene[antRes[[2]]$abdomen=="queen"]]="queen"
ogg11$abdBee[ogg11$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="worker"]]="worker"
ogg11$abdBee[ogg11$gene_Amel %in% beeRes[[2]]$Gene[beeRes[[2]]$abdomen=="queen"]]="queen"
ogg11$BeeQueen = ogg11$BeeWorker = ogg11$AntQueen = ogg11$AntWorker = ogg11$AntCaste = ogg11$BeeCaste = ogg11$conCaste = ogg11$conDev = ogg11$conQueen = ogg11$conWorker = 0
ogg11$BeeQueen[ogg11$abdBee=="queen"]=1
ogg11$BeeWorker[ogg11$abdBee=="worker"]=1
ogg11$AntQueen[ogg11$abdAnt=="queen"]=1
ogg11$AntWorker[ogg11$abdAnt=="worker"]=1
ogg11$AntCaste[ogg11$abdAnt!="nonDE"]=1
ogg11$BeeCaste[ogg11$abdBee!="nonDE"]=1

ogg11$conQueen[ogg11$AntQueen + ogg11$BeeQueen == 2] = 1
ogg11$conWorker[ogg11$AntWorker + ogg11$BeeWorker == 2] = 1
ogg11$conCaste[ogg11$AntCaste + ogg11$BeeCaste == 2] = 1
ogg11$conDev[ogg11$DevAnt + ogg11$DevBee == 2] = 1

x <- vennCounts(ogg11[,c("conDev","conWorker")])

#Figure 2b. 
png("~/GitHub/devnetwork/figures/DevCaste.png")
vennDiagram(x,names = c("development","caste"))
dev.off()

####
devOverlap <- function(devGene,DEres){
  d = data.frame(Gene = rownames(DEres),Dev=0,DE=0)
  d$Dev[d$Gene %in% devGene]=1
  d$DE[DEres$FDR < 0.05]=1
  return(d)
}

devFisher <- function(devGene,tests){
  ACD <- lapply(tests,function(x) devOverlap(devGene,x))
  ACDvenn <- lapply(ACD,function(x) vennCounts(x[,c(2,3)]))
  AvennChi <- lapply(ACDvenn,function(x) fisher.test(rbind(c(x[4,3],x[3,3]),c(x[2,3],x[1,3]))))
  return(list(ACDvenn,AvennChi))
}

consDevOverlap <- function(d,col1,col2,name1,name2){
  d = d[,colnames(d) %in% c("DevAnt","DevBee",col1,col2)]
  d$consDev=0
  d$consDev[d$DevAnt+d$DevBee==2]=1
  d$consCaste = d$consQ = d$consW = d$cons_reverse = 0
  d$consCaste[d[,col1]!="nonDE" & d[,col2] != "nonDE"]=1
  d$consQ[d[,col1]==name1 & d[,col2] == name1]=1
  d$consW[d[,col1]==name2 & d[,col2] == name2]=1
  d$cons_reverse[(d[,col1]==name2 & d[,col2] == name1) | (d[,col1]==name1 & d[,col2] == name2)]=1
  return(d)
}

ACD = devFisher(AntDev,antTests_oneLarv[c(3:5)])
BCD = devFisher(BeeDev,beeTests_oneLarv[c(3:5)])
ASD = devFisher(AntDev,antSocial)
BSD = devFisher(BeeDev,beeSocial)

plotOverlap <- function(dL){
  
}

d = ogg11
d = d[,colnames(d) %in% c("OGGend","gene_Mphar","gene_Amel","DevAnt","DevBee")]
d = merge(d,antRes[[2]],by.x="gene_Mphar",by.y="Gene")
d = merge(d,beeRes[[2]],by.x="gene_Amel",by.y="Gene")
CasteDev <- lapply(c("head","thorax","abdomen"),function(x) consDevOverlap(d,paste(x,"x",sep="."),paste(x,"y",sep="."),"queen","worker"))
CasteDev <- lapply(CasteDev,function(x) x[,c(5:9)])
CasteDev_overlap <- lapply(CasteDev,function(x){
  lapply(c(2:5), function(j) vennCounts(x[,c(1,j)]))
})
CasteDev_fisher <- lapply(unlist(CasteDev_overlap,recursive=FALSE),function(x) fisher.test(rbind(c(x[4,3],x[3,3]),c(x[2,3],x[1,3]))))

d = ogg11
d = d[,colnames(d) %in% c("OGGend","gene_Mphar","gene_Amel","DevAnt","DevBee")]
d = merge(d,antSocRes[[2]],by.x="gene_Mphar",by.y="Gene")
d = merge(d,beeSocRes[[2]],by.x="gene_Amel",by.y="Gene")
CasteDev <- lapply(c("head","mesosoma","gaster"),function(x) consDevOverlap(d,paste(x,"x",sep="."),paste(x,"y",sep="."),"nurse","forager"))
CasteDev <- lapply(CasteDev,function(x) x[,c(5:9)])
CasteDev_overlap <- lapply(CasteDev,function(x){
  lapply(c(2:5), function(j) vennCounts(x[,c(1,j)]))
})
CasteDev_fisher <- lapply(unlist(CasteDev_overlap,recursive=FALSE),function(x) fisher.test(rbind(c(x[4,3],x[3,3]),c(x[2,3],x[1,3]))))


oggCaste <- function(a,b){
  Aq = ogg11$OGG[ogg11$gene_Mphar %in% a[[2]]]
  Anq = ogg11$OGG[ogg11$gene_Mphar %in% a[[3]]]
  Bq = ogg11$OGG[ogg11$gene_Amel %in% b[[2]]]
  Bnq = ogg11$OGG[ogg11$gene_Amel %in% b[[3]]]
  qC = Aq[Aq %in% Bq]
  nqC = Anq[Anq %in% Bnq]
  bee_cons_q = ogg11$gene_Amel[ogg11$OGG %in% qC]
  bee_cons_nq = ogg11$gene_Amel[ogg11$OGG %in% nqC]
  ant_cons_q = ogg11$gene_Mphar[ogg11$OGG %in% qC]
  ant_cons_nq = ogg11$gene_Mphar[ogg11$OGG %in% nqC]
  return(list(list(ant_cons_q,bee_cons_q),list(ant_cons_nq,bee_cons_nq)))
}

antCaste = ogg11$gene_Mphar[ogg11$abdBee!="nonDE" & ogg11$abdAnt!="nonDE"]
beeCaste = ogg11$gene_Amel[ogg11$abdBee!="nonDE" & ogg11$abdAnt!="nonDE"]

ant_sexDE_filt = lapply(ant_sexDE,function(x) x[rownames(x) %in% antCaste,])
AsexRes <- lapply(ant_sexDE_filt,extractBias)
antCasteDE_filt = lapply(antTests_oneLarv[c(3:5)],function(x) x[rownames(x) %in% antCaste,] )
AcasteRes <- lapply(antCasteDE_filt,extractBias)
AFC <- lapply(c(1:3),function(x) FCplot_filt(AsexRes[[x]],AcasteRes[[x]]))

bee_sexDE_filt = lapply(bee_sexDE,function(x) x[rownames(x) %in% beeCaste,])
BsexRes <- lapply(bee_sexDE_filt,extractBias)
beeCasteDE_filt = lapply(beeTests_oneLarv[c(3:5)],function(x) x[rownames(x) %in% beeCaste,] )
BcasteRes <- lapply(beeCasteDE_filt,extractBias)
BFC <- lapply(c(1:3),function(x) FCplot_filt(BsexRes[[x]],BcasteRes[[x]]))

png("~/GitHub/devnetwork/figures/Caste_sexFC.png",height=2000,width=4000,res=300)
grid.arrange(
AFC[[3]][[1]]+geom_smooth(aes(color=DE),method="lm",se=FALSE)+xlab("sex bias")+ylab("caste bias")+ggtitle("ant"),
BFC[[3]][[1]]+geom_smooth(aes(color=DE),method="lm",se=FALSE)+xlab("sex bias")+ylab("caste bias")+ggtitle("bee"),
nrow=1)
dev.off()

AsexG = ant_sexDE[[3]]
AsexG$Gene = rownames(AsexG)
AcasteG = antTests_oneLarv[[5]]
AcasteG$Gene = rownames(AcasteG)
BsexG = bee_sexDE[[3]]
BsexG$Gene = rownames(BsexG)
BcasteG = beeTests_oneLarv[[5]]
BcasteG$Gene = rownames(BcasteG)
ogg11 = merge(ogg11,AsexG[,c(1,6)],by.x="gene_Mphar",by.y="Gene")
ogg11 = merge(ogg11,AcasteG[,c(1,6)],by.x="gene_Mphar",by.y="Gene")
ogg11 = merge(ogg11,BsexG[,c(1,6)],by.x="gene_Amel",by.y="Gene")
ogg11 = merge(ogg11,BcasteG[,c(1,6)],by.x="gene_Amel",by.y="Gene")
colnames(ogg11)[c(6:9)] = c("sexAnt","casteAnt","sexBee","casteBee")
ogg11$casteDE = "nonDE"
ogg11$casteDE[ogg11$abdBee!="nonDE" & ogg11$abdAnt != "nonDE"]="flipped"
ogg11$casteDE[ogg11$abdBee=="queen" & ogg11$abdAnt == "queen"]="queen_both"
ogg11$casteDE[ogg11$abdBee=="worker" & ogg11$abdAnt == "worker"]="worker_both"


p1 <- ggplot(ogg11[ogg11$casteDE!="nonDE",],aes(x = -sexAnt,y=-casteAnt,color=casteDE))+
  geom_point(alpha=0.5)+main_theme+
  theme(legend.position = c(0.8,0.2))+
  xlab("logFC (ant queen/ant male)")+
  ylab("logFC (ant queen/ant worker)")

p2 <- ggplot(ogg11[ogg11$casteDE!="nonDE",],aes(x = -sexBee,y=-casteBee,color=casteDE))+
  geom_point(alpha=0.5)+main_theme+
  theme(legend.position = c(0.8,0.2))+
  xlab("logFC (bee queen/bee male)")+
  ylab("logFC (bee queen/bee worker)")

p3 <- ggplot(ogg11[ogg11$casteDE!="nonDE",],aes(x = -sexBee,y=-sexAnt,color=casteDE))+
  geom_point(alpha=0.5)+main_theme+
  theme(legend.position = c(0.8,0.2))+
  xlab("logFC (bee queen/bee male)")+
  ylab("logFC (ant queen/ant male)")

png("~/GitHub/devnetwork/figures/Caste_Sex.png",width=5000,height=2000,res=300)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

ant_OGG = merge(ogg11[,c(1:5)],AsexG[,c(1,6)],by.x="gene_Mphar",by.y="Gene",all.y=T)
ant_OGG = merge(ant_OGG,AcasteG[,c(1,6)],by.x="gene_Mphar",by.y="Gene",all.y=T)
ant_OGG = ant_OGG[ant_OGG$gene_Mphar %in% rownames(antTests_oneLarv[[5]])[antTests_oneLarv[[5]]$FDR < 0.05],]
colnames(ant_OGG)[c(6,7)] = c("antSex","antCaste")
ant_OGG$casteDE = "LS_noOGG"
ant_OGG$casteDE[!is.na(ant_OGG$abdBee)] = "LS_OGG"
ant_OGG$casteDE[ant_OGG$abdBee!="nonDE" & !is.na(ant_OGG$abdBee)] = "consDE"

fac = "consDE"
cor.test(ant_OGG$antSex[ant_OGG$casteDE==fac],ant_OGG$antCaste[ant_OGG$casteDE==fac])
ggplot(ant_OGG,aes(x = antSex,y=antCaste,color=casteDE))+
  geom_point(alpha=0.5)


AvmRes <- lapply(ant_VM,extractBias)
AsocRes <- lapply(antSocial,extractBias)
AsexRes <- lapply(ant_sexDE,extractBias)
AcasteRes <- lapply(antTests_oneLarv[c(3:5)],function(x) extractBias(x))

AFC <- lapply(c(1:3),function(x) FCplot(AsexRes[[x]],AcasteRes[[x]]))
AFC_vm <- lapply(c(1:3),function(x) FCplot(AvmRes[[x]],AcasteRes[[x]]))
AFC_NF_caste <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AcasteRes[[x]]))
AFC_NF_sex <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AsexRes[[x]]))
AFC_NF_vm <- lapply(c(1:3),function(x) FCplot(AsocRes[[x]],AvmRes[[x]]))

BsexRes <- lapply(bee_sexDE,extractBias)
BcasteRes <- lapply(beeTests_oneLarv[c(3:5)],function(x) extractBias(x))
BvmRes <- lapply(bee_VM,extractBias)
BsocRes <- lapply(beeSocial,extractBias)
BFC <- lapply(c(1:3),function(x) FCplot(BsexRes[[x]],BcasteRes[[x]]))
BFC_vm <- lapply(c(1:3),function(x) FCplot(BvmRes[[x]],BcasteRes[[x]]))
BFC_NF_caste <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BcasteRes[[x]]))
BFC_NF_sex <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BsexRes[[x]]))
BFC_NF_vm <- lapply(c(1:3),function(x) FCplot(BsocRes[[x]],BvmRes[[x]]))

AFCplot <- lapply(AFC,function(x) x[[1]]+ylab("queen/worker log2 FC")+xlab("queen/male log2 FC"))
legend <- g_legend(AFCplot[[1]]+theme(legend.position="right"))

png("~/GitHub/devnetwork/figures/scLFCant.png",height=2000,width=2000,res=300)
grid.arrange(AFCplot[[1]]+theme(legend.position="hidden")+ggtitle("ant head"),
             AFCplot[[2]]+theme(legend.position="hidden")+ggtitle("ant thorax"),
             AFCplot[[3]]+theme(legend.position="hidden")+ggtitle("ant abdomen"),legend,ncol=2)
dev.off()

AFCplot <- lapply(AFC_NF_sex,function(x) x[[1]]+ylab("nurse/forager log2 FC")+xlab("queen/male log2 FC"))
png("~/GitHub/devnetwork/figures/socLFCant.png",height=2000,width=4000,res=300)
grid.arrange(AFCplot[[1]]+theme(legend.position="none")+ggtitle("head"),
             AFCplot[[2]]+theme(legend.position="none")+ggtitle("thorax"),
             AFCplot[[3]]+theme(legend.position="none")+ggtitle("abdomen"),ncol=3)
dev.off()

AFCplot <- lapply(BFC_NF_sex,function(x) x[[1]]+ylab("nurse/forager log2 FC")+xlab("queen/male log2 FC"))
png("~/GitHub/devnetwork/figures/socLFCbee.png",height=2000,width=4000,res=300)
grid.arrange(AFCplot[[1]]+theme(legend.position="none")+ggtitle("head"),
             AFCplot[[2]]+theme(legend.position="none")+ggtitle("thorax"),
             AFCplot[[3]]+theme(legend.position="none")+ggtitle("abdomen"),ncol=3)
dev.off()

BFCplot <- lapply(BFC,function(x) x[[1]]+ylab("queen/worker log2 FC")+xlab("queen/male log2 FC"))
legend <- g_legend(BFCplot[[1]]+theme(legend.position="right"))

png("~/GitHub/devnetwork/figures/scLFCbee.png",height=2000,width=4000,res=300)
grid.arrange(BFCplot[[1]]+theme(legend.position="none")+ggtitle("head"),
             BFCplot[[2]]+theme(legend.position="none")+ggtitle("thorax"),
             BFCplot[[3]]+theme(legend.position="none")+ggtitle("abdomen"),ncol=3)
dev.off()

consCaste <- mapply(oggCaste,AcasteRes,BcasteRes,SIMPLIFY=FALSE)
consSex <- mapply(oggCaste,AsexRes,BsexRes,SIMPLIFY=FALSE)
consVM <- mapply(oggCaste,AvmRes,BvmRes,SIMPLIFY=FALSE)

vennConserved <- function(caste,sex,vm,queen){
  o = ogg11
  o$Caste = o$Sex = o$VM = 0
  o$Caste[o$gene_Mphar %in% caste[[queen]][[1]]] = 1
  o$Sex[o$gene_Mphar %in% sex[[queen]][[1]]] = 1
  o$VM[o$gene_Mphar %in% vm[[queen]][[1]]] = 1
  x <- vennCounts(o[,c("Caste","Sex","VM")])
  return(list(x,o))
}

Cq <- lapply(c(1:3),function(x) vennConserved(consCaste[[x]],consSex[[x]],consVM[[x]],1))
Cw <- lapply(c(1:3),function(x) vennConserved(consCaste[[x]],consSex[[x]],consVM[[x]],2))

AntSC = merge(AcasteRes[[3]][[1]],AsexRes[[3]][[1]],by = "Gene")
AntSC$conCaste = 0
AntSC$conCaste[AntSC$Gene %in% c(as.character(consCaste[[3]][[1]][[1]]),as.character(consCaste[[3]][[2]][[1]]))] = 1
ggplot(AntSC,aes(x=-FC.x,y=-FC.y,color=as.factor(conCaste)))+
  geom_point()+geom_smooth()

BeeSC = merge(BcasteRes[[3]][[1]],BsexRes[[3]][[1]],by = "Gene")
BeeSC$conCaste = 0
BeeSC$conCaste[BeeSC$Gene %in% c(as.character(consCaste[[3]][[1]][[2]]),as.character(consCaste[[3]][[2]][[2]]))] = 1
ggplot(BeeSC,aes(x=-FC.x,y=-FC.y,color=as.factor(conCaste)))+
  geom_point()+geom_smooth()


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

cbScatter <- function(CB,col1,col2){
  p <- ggplot(CB,aes_q(x=as.name(col1),
              y=as.name(col2)))+
    geom_point(alpha = 0.5,size = 0.7)+
    geom_smooth(se=FALSE,method="lm")+
    main_theme
  test = cor.test(CB[,col1],CB[,col2],method="spearman")
  return(list(p,test))
}

antCB = euclDist(antRes_allstage[[1]])
beeCB = euclDist(beeRes_allstage[[1]])

antCBps = merge(antCB,antPS,by="Gene")
cor.test(antCBps$cb_noAdult,antCBps$ps)
ggplot(antCBps,aes(x=psName,y=cb_noAdult))+
  geom_boxplot()


antCBS = data.frame(Gene = antSocRes[[1]]$Gene,Soc_bias=apply(antSocRes[[1]][,-c(1)],1,function(x) sqrt(sum(x^2))/length(x)))
beeCBS = data.frame(Gene = beeSocRes[[1]]$Gene,Soc_bias=apply(beeSocRes[[1]][,-c(1)],1,function(x) sqrt(sum(x^2))/length(x)))
ant_bias = merge(antCB,antCBS,by="Gene")
bee_bias = merge(beeCB,beeCBS,by="Gene")

png("~/GitHub/devnetwork/figures/cb_cor.png",width=4000,height=5000,res=300)
grid.arrange(cbScatter(ant_bias,"cb_larva","cb_adult")[[1]]+ggtitle("ant"),
             cbScatter(bee_bias,"cb_larva","cb_adult")[[1]]+ggtitle("bee"),
             cbScatter(ant_bias,"cb_larva","Soc_bias")[[1]]+ggtitle("ant"),
             cbScatter(bee_bias,"cb_larva","Soc_bias")[[1]]+ggtitle("bee"),
             cbScatter(ant_bias,"cb_adult","Soc_bias")[[1]]+ggtitle("ant"),
             cbScatter(bee_bias,"cb_adult","Soc_bias")[[1]]+ggtitle("bee"),nrow=3)
dev.off()

ABcb = merge(ant_bias,AB11,by.x = "Gene",by.y="gene_Mphar")
ABcb = merge(ABcb,bee_bias,by.x="gene_Amel",by.y="Gene")
ABcb = merge(ABcb,melt(antDE,id.vars="Gene"),by="Gene")
ABcb = merge(ABcb,melt(beeDE,id.vars="Gene"),by.x=c("variable","gene_Amel"),by.y=c("variable","Gene"))

ABcb$cons = "nonDE"
ABcb$cons[ABcb$value.x!=ABcb$value.y] = "inconsistent"
ABcb$cons[ABcb$value.x=="worker" & ABcb$value.y=="worker"] = "worker"
ABcb$cons[ABcb$value.x=="queen" & ABcb$value.y=="queen"] = "queen"
ABcb$casteDE = ABcb$cons
ABcb$casteDE[ABcb$value.x=="worker" & ABcb$value.y=="queen"] = "AW-BQ"
ABcb$casteDE[ABcb$value.x=="worker" & ABcb$value.y=="nonDE"] = "AW-NDE"
ABcb$casteDE[ABcb$value.x=="queen" & ABcb$value.y=="worker"] = "AQ-BW"
ABcb$casteDE[ABcb$value.x=="queen" & ABcb$value.y=="nonDE"] = "AQ-NDE"
ABcb$casteDE[ABcb$value.x=="nonDE" & ABcb$value.y=="queen"] = "NDE-BQ"
ABcb$casteDE[ABcb$value.x=="nonDE" & ABcb$value.y=="worker"] = "NDE-BW"

png("~/GitHub/devnetwork/figures/cb_cor_species.png",width=6000,height=2000,res=300)
grid.arrange(cbScatter(ABcb,"cb_larva.x","cb_larva.y")[[1]],
             cbScatter(ABcb,"cb_adult.x","cb_adult.y")[[1]],
             cbScatter(ABcb,"Soc_bias.x","Soc_bias.y")[[1]],nrow=1)
dev.off()

allCB = merge(ant_bias,AB11,by.x = "Gene",by.y="gene_Mphar")
allCB$species = "Mphar"
allCB2 = merge(bee_bias,AB11,by.x = "Gene",by.y="gene_Amel")
allCB2$species="Amel"
allCB3 = rbind(allCB[,-c(10)],allCB2[-c(10)])
allCB3 = merge(allCB3,ABcb[,c("OGGend","cons","variable","casteDE")],by = "OGGend")
allCB3$cons = factor(allCB3$cons,levels = c("nonDE","inconsistent","queen","worker"))

png("~/GitHub/devnetwork/figures/DEtype_castebias.png",width=2000,height=2000,res=300)
ggplot(allCB3[allCB3$variable=="abdomen",],aes(y=cb_larva,x=species,fill = cons))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  xlab("DE type")+
  ylab("larval caste bias")+
  ylim(0,1.25)+
  main_theme
dev.off()

png("~/GitHub/devnetwork/figures/DEtype_castebias.png",width=2000,height=2000,res=300)
ggplot(allCB3[allCB3$variable=="abdomen",],aes(y=cb_larva,x=species,fill = casteDE))+
  geom_boxplot(notch=TRUE,outlier.shape=NA)+
  xlab("DE type")+
  ylab("larval caste bias")+
  ylim(0,1.25)+
  main_theme
dev.off()

ant_allFC = merge(antRes_allstage[[1]],antSocRes[[1]],by="Gene")

bee_allFC = merge(beeRes_allstage[[1]],beeSocRes[[1]],by="Gene")

png("~/GitHub/devnetwork/figures/cb_cor.png",width=5000,height=4000,res=300)
grid.arrange(cbScatter(ant_allFC,"8_gaster","gaster")[[1]]+ggtitle("ant abd")+ylab("for/nurse logFC")+xlab("W/Q logFC"),
             cbScatter(ant_allFC,"8_mesosoma","mesosoma")[[1]]+ggtitle("ant thorax")+ylab("for/nurse logFC")+xlab("W/Q logFC"),
             cbScatter(ant_allFC,"8_head","head")[[1]]+ggtitle("ant head")+ylab("for/nurse logFC")+xlab("W/Q logFC"),
             cbScatter(bee_allFC,"8_gaster","gaster")[[1]]+ggtitle("bee abd")+ylab("for/nurse logFC")+xlab("W/Q logFC"),
             cbScatter(bee_allFC,"8_mesosoma","mesosoma")[[1]]+ggtitle("bee thorax")+ylab("for/nurse logFC")+xlab("W/Q logFC"),
             cbScatter(bee_allFC,"8_head","head")[[1]]+ggtitle("bee head")+ylab("for/nurse logFC")+xlab("W/Q logFC"),nrow=2)
dev.off()

#Load ps 
Amel_acu <- read.table("~/GitHub/devnetwork/phylo_results/Amel_acu")
Amel_blast <- read.table("~/Gi")





x <- vennCounts(ogg11[,c("DevBee","DevAnt")])

#Figure 2a
png("~/GitHub/devnetwork/figures/DevelomentalOverlap.png")
vennDiagram(x,names = c("Apis","Mphar"))
dev.off()


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

Dmel_oggDat = merge(DmelDat,AB111,by="Gene")
Dogg = merge(Dmel_oggDat,ogg11,by="gene_Mphar")

Dabd = merge(ogg11,Dmel_oggDat,by="OGGend")
Dabd$casteDE[Dabd$abdBee=="queen" & Dabd$abdAnt=="worker"] = "BQ-AW"
Dabd$casteDE[Dabd$abdBee=="worker" & Dabd$abdAnt=="queen"] = "BW-AQ"
Dabd$casteDE[Dabd$abdBee=="worker" & Dabd$abdAnt=="nonDE"] = "BW-NDE"
Dabd$casteDE[Dabd$abdBee=="queen" & Dabd$abdAnt=="nonDE"] = "BQ-NDE"
Dabd$casteDE[Dabd$abdBee=="nonDE" & Dabd$abdAnt=="worker"] = "AW-NDE"
Dabd$casteDE[Dabd$abdBee=="nonDE" & Dabd$abdAnt=="queen"] = "AQ-NDE"

png("~/GitHub/devnetwork/figures/DE_dmelCasteLogFC.png",width=2000,height=2000,res=300)
ggplot(Dabd,aes(x = casteDE,y=-logFC))+
  geom_violin()+
  ylab("logFC (dmel female/dmel male)")+
  geom_boxplot(notch=TRUE,width=0.2,outlier.shape = NA)+
  main_theme+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        plot.margin = unit(c(0.5,1.4,0.5,0.5),"cm"))
dev.off()


AcasteDmel <- lapply(antTests_oneLarv[c(3:5)],function(x){
  x$Gene = rownames(x)
  return(merge(x,Dmel_oggDat,by.x="Gene",by.y="gene_Mphar"))
})
AdPlot <- lapply(AcasteDmel,DmelSexFC)

Adplot <- lapply(AdPlot,function(x) x[[1]]+ylab("ant queen/worker log2 FC")+xlab("Dmel female/male log2 FC"))

BcasteDmel <- lapply(beeTests_oneLarv[c(3:5)],function(x){
  x$Gene = rownames(x)
  return(merge(x,Dmel_oggDat,by.x="Gene",by.y="gene_Amel"))
})
BdPlot <- lapply(BcasteDmel,DmelSexFC)

Bdplot <- lapply(BdPlot,function(x) x[[1]]+ylab("bee queen/worker log2 FC")+xlab("Dmel female/male log2 FC"))

png("~/GitHub/devnetwork/figures/scLFCant_dmel.png",height=2000,width=4000,res=300)
grid.arrange(Adplot[[1]]+theme(legend.position="none")+ggtitle("head"),
             Adplot[[2]]+theme(legend.position="none")+ggtitle("thorax"),
             Adplot[[3]]+theme(legend.position="none")+ggtitle("abdomen"),ncol=3)
dev.off()

png("~/GitHub/devnetwork/figures/scLFCbee_dmel.png",height=2000,width=4000,res=300)
grid.arrange(Bdplot[[1]]+theme(legend.position="none")+ggtitle("head"),
             Bdplot[[2]]+theme(legend.position="none")+ggtitle("thorax"),
             Bdplot[[3]]+theme(legend.position="none")+ggtitle("abdomen"),ncol=3)
dev.off()


DmelSexFC <- function(data){
  data$DE = "nonDE"
  data$DE[data$logFC.x < 0 & data$logFC.y < 0 & data$FDR.x < 0.05 & data$FDR.y < 0.05] = "female-upregulated"
  data$DE[data$logFC.x > 0 & data$logFC.y > 0 & data$FDR.x < 0.05 & data$FDR.y < 0.05] = "female-downregulated"
  p <- ggplot(data,aes(x = -logFC.y,y = -logFC.x))+ #Queen-up will be positive
    geom_point(aes(color = DE),alpha=0.5)+
    geom_smooth(method="lm",se=FALSE,color="black")+
    scale_color_manual(values = c("blue","red","grey60"))+
    main_theme
  return(list(p,cor.test(data$logFC.x,data$logFC.y,method = "spearman"),table(data$DE)))
}

dmelF_B <- Dmel_oggDat$gene_Amel[Dmel_oggDat$FDR < 0.05 & Dmel_oggDat$logFC < 0]
dmelF_A <- Dmel_oggDat$gene_Amel[Dmel_oggDat$FDR < 0.05 & Dmel_oggDat$logFC < 0]

Acaste <- c(BcasteRes[[3]][[3]])
Acaste = Acaste[Acaste %in% Dmel_oggDat$gene_Amel]
AoverQ <- sum(dmelF_A %in% Acaste)

t <- rbind(c(AoverQ,length(dmelF_A) - AoverQ),
           c(length(Acaste) - AoverQ,
             nrow(Dmel_oggDat) - length(Acaste) - length(dmelF_A) + AoverQ))

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

antT = antT[rowSums(antT) > 0,]
Agenes = data.frame(Gene=rownames(antT),index = seq(0,(nrow(antT) - 1)))
beeT = beeT[rowSums(beeT) > 0,]
Bgenes = data.frame(Gene=rownames(beeT),index = seq(0,(nrow(beeT) - 1)))
AB1 = merge(AB11,Agenes,by.x="gene_Mphar",by.y="Gene")
AB1 = merge(AB1,Bgenes,by.x="gene_Amel",by.y="Gene")
AB1$a1 = 1
AB1$a2 = 2
AB1$a3 = 1
write.table(AB1[,c(6:8,4,5)],file="~/GitHub/devnetwork/data/orthologymap_integers.txt",col.names=FALSE,row.names = FALSE,sep='\t')


orthologs <- merge(Aps[!is.na(Aps$OGGacu),c("Gene","OGGacu")],Bps[!is.na(Aps$OGGacu),c("Gene","OGGacu")],by="OGGacu")
colnames(orthologs) = c("OGGacu","gene_Mphar","gene_Amel")
#Get ortholog weighting values for use 

AmelIndex = data.frame(Gene = rownames(beeT_filt),index = seq(0,nrow(beeT_filt) - 1))
MpharIndex = data.frame(Gene = rownames(antT_filt),index = seq(0,nrow(antT_filt) - 1))
orthologs = merge(orthologs,MpharIndex,by.x="gene_Mphar",by.y="Gene")
orthologs = merge(orthologs,AmelIndex,by.x="gene_Amel",by.y="Gene")
t = table(droplevels(orthologs$OGGacu))
weight = data.frame(OGGacu = names(t), mag = 1/as.numeric(as.character(t)))
orthologs = merge(orthologs,weight,by="OGGacu")
OGGmap = orthologs[,c(4,5,1,6)]
write.table(OGGmap,file="~/Data/devnetwork/OGGmap.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)

########Plaid clustering
library(biclust)

plaid <- function(data,boots){
  res = list()
  for (i in 1:boots){
    print(i)
    res[[i]] = biclust(as.matrix(log(data+(data^2+1)^0.5)),method=BCPlaid(),max.layers=25,iter.startup=10,iter.layer=15,verbose=FALSE)
  }
  return(res) 
}

antPl = plaid(antT,100)
beePl = plaid(beeT,100)

#using an inverse hyperbolic sine tranformation
b <- biclust(as.matrix(log(beeT+(beeT^2+1)^0.5)),method=BCPlaid(),max.layers=25,iter.startup=10,iter.layer=15)
#Bicluster 7 is for adult queen gaster and eggs
bic <- b@RowxNumber[,7]
genes = rownames(beeT)[bic]
testB = beeRes[[2]]
testB$bic = 0
testB$bic[testB$Gene %in% genes]=1

ACUogg$bicBee = 0
ACUogg$bicBee[ACUogg$gene_Amel %in% genes] = 1

ggplot(testB,aes(x=as.factor(bic),y=abdomen))+
  geom_boxplot()

b <- biclust(as.matrix(log(antT+(antT^2+1)^0.5)),method=BCPlaid(),max.layers=25,iter.startup=10,iter.layer=15)
#Bicluster 7 is for adult queen gaster and eggs
bic <- b@RowxNumber[,9]
genes = rownames(antT)[bic]
testB = antRes[[2]]
testB$bic = 0
testB$bic[testB$Gene %in% genes]=1

ggplot(testB,aes(x=as.factor(bic),y=abdomen))+
  geom_boxplot(notch=TRUE)



