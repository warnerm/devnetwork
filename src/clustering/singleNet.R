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
  write.csv(tCor[1:1000,1:1000],paste("~/Data/devnetwork/",species,"TESTpCor.csv",sep=""),row.names=FALSE)
  write.csv(tCor,paste("~/Data/devnetwork/",species,"pCor.csv",sep=""),row.names=FALSE)
  return(0)
}

makeCor("bees")
makeCor("ants")

###Prepare data for OrthoClust
prepOrthoClust <- function(species,d){
  df <- read.csv(paste("~/Data/devnetwork/",species,"TESTpCor.csv",sep=""))
  df = abs(df)
  ranks <- apply(df,1,order,decreasing=TRUE)
  input_df <- ranks[2:(d+1),]
  input <- matrix(ncol=2)
  for (i in 1:nrow(input_df)){
    for (j in 1:ncol(input_df)){
      input <- rbind(input,c(j,input_df[i,j])) #Column is the actual value of the gene, row is 
    }
  }
}
