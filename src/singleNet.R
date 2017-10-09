#Will use IndivNet.py to determine single-species networks
#This script takes in tpm matrices and outputs pearson correlation matrices
makeCor <- function(species){
  tpm <- read.table(paste("~/GitHub/devnetwork/data/",species,".tpm.txt",sep=""),header=TRUE)
  tpm <- tpm[!grepl("ERCC",tpm$id),]
  rownames(tpm) = tpm$id
  tpm = tpm[,-c(1)]
  tpm = tpm[rowSums(tpm) > 0,]
  tCor = cor(t(tpm))
  tCor = round(tCor,2)
  write.csv(tCor,paste("~/GitHub/devnetwork/data/data.processed/",species,"pCor.csv",sep=""))
  return(0)
}

makeCor("bees")
makeCor("ants")