load("../data/DEtests.RData")
beeMK <- read.csv("../out/Mkinput.csv",header = F)
beeSub <- read.csv("../out/substitutions.csv")
beeMK = cbind(beeMK,rep(1,nrow(beeMK)))

antEvol <- read.csv("../data/MpharAnn.csv")
antSub = antEvol[!is.na(antEvol$Fixed.Non.Synonymous),c(1,13,17,15,17,12,16,14,16)]

##22 chromosomes sampled
antMK = cbind(antSub[,-c(1)],rep(22,nrow(antSub)),rep(1,nrow(antSub)))

addClass <- function(DEdat,sub,mk,species){
  mk$V10 = rep(1,nrow(mk))
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="queen"]] = 2
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="worker"]] = 3
  write.table(mk,file = paste("../MK_alpha_input/",names(DEdat)[2],".",species,".csv",sep = ""),col.names = FALSE,row.names = FALSE,sep=",")
}

sapply(c(2:ncol(beeRes_allstage[[2]])),function(x) addClass(beeRes_allstage[[2]][,c(1,x)],beeSub,beeMK,"bee"))
sapply(c(2:ncol(antRes_allstage[[2]])),function(x) addClass(antRes_allstage[[2]][,c(1,x)],antSub,antMK,"ant"))