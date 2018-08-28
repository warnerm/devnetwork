library(rslurm)

plaid <- function(run){
  library(biclust)
  print(run)
  a = biclust(as.matrix(log(data+(data^2+1)^0.5)),method=BCPlaid(),max.layers=25,iter.startup=10,iter.layer=15,verbose=FALSE)
  return(a)
}

#Take text file and get data frame we want for analysis
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}

beeT <- read.table("../data/bees.tpm.txt",header=TRUE)
antT <- read.table("../data/ants.tpm.txt",header=TRUE)
beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,]
beeT = beeT[rowSums(beeT) > 0,]




data = antT
runs = data.frame(run = seq(1,100,by=1))
sjob <- slurm_apply(plaid, runs, jobname = 'plaid',
                    nodes = 4, cpus_per_node = 3, add_objects="data",submit = TRUE)
antPl <- get_slurm_out(sjob,outtype='raw') #get output as lists


data = beeT
runs = data.frame(run = seq(1,100,by=1))
sjob <- slurm_apply(plaid, runs, jobname = 'plaid',
                    nodes = 4, cpus_per_node = 3, add_objects="data",submit = TRUE)
beePl <- get_slurm_out(sjob,outtype='raw') #get output as lists


save(antPl,beePl,file="PlaidResults.RData")