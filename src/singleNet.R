#Will use IndivNet.py to determine single-species networks
#This script takes in tpm matrices and outputs pearson correlation matrices
makeCor <- function(species){
  tpm <- read.table(paste("~/GitHub/devnetwork/data/",species,".tpm.txt",sep=""),header=TRUE)
  tpm <- tpm[!grepl("ERCC",tpm$id),]
  rownames(tpm) = tpm$id
  tpm = tpm[,-c(1)]
  tpm = tpm[rowSums(tpm) > 0,]
  tCor = cor(t(tpm))
  tCor = round(tCor,3)
  write.csv(tCor[1:100,1:100],paste("~/Data/devnetwork/",species,"TESTpCor.csv",sep=""),row.names=FALSE)
  write.csv(tCor,paste("~/Data/devnetwork/",species,"pCor.csv",sep=""),row.names=FALSE)
  return(0)
}

makeCor("bees")
makeCor("ants")