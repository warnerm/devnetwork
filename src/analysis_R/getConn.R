setwd("~/GitHub/devnetwork/")
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}
beeT <- read.table("data/bees.tpm.txt",header=TRUE)
antT <- read.table("data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]


#Calculate connectivity 
antCor = cor(t(antT[rowSums(antT) > 0,]))^6
beeCor = cor(t(beeT[rowSums(beeT) > 0,]))^6
antConn = data.frame(Gene = rownames(antCor),kTotal = rowSums(antCor))
beeConn = data.frame(Gene = rownames(beeCor), kTotal = rowSums(beeCor))

write.csv(antConn,"results/antConnectivity.csv")
write.csv(beeConn,"results/beeConnectivity.csv")