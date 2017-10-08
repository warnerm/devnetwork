library(edgeR)
load("initialvariables.RData")
##bootstrapping the whole process
bootOverlap <- function(antF,beeF,antC,beeC,test){
  bee <- list()
  ant <- list()
  TwoSpec <- list()
  
  for (tissue in tissues){
    if (test == "caste"){
      bee[[tissue]]=tissueCaste(beeF,beeC,tissue,scramble=TRUE)
      ant[[tissue]]=tissueCaste(antF,antC,tissue,scramble=TRUE)
    } else {
      bee[[tissue]]=tissueSocial(beeF,beeC,tissue,scramble=TRUE)
      ant[[tissue]]=tissueSocial(antF,antC,tissue,scramble=TRUE)
    }
    beeOGG <- ogg2$OGG[ogg2$gene_Amel %in% bee[[tissue]]]
    antOGG <- ogg2$OGG[ogg2$gene_Mphar %in% ant[[tissue]]]
    TwoSpec[[tissue]] <- beeOGG[beeOGG %in% antOGG]
  }
  
  AllTest <- list(bee,ant,TwoSpec)
  return(calcOverlap(AllTest))
}

#Calculate overlap of tests
calcOverlap <- function(AllTest){
  overlap = matrix(nrow=1,ncol=9)
  for (i in 1:3){
    for (j in 1:3){
      overlap[1,(i-1)*3+j] = sum(AllTest[[i]][[tissues[j]]] %in% DevList[[i]])
    }
  }
  colnames(overlap) = apply(expand.grid(tissues,c("bee","ant","overlap")), 1, paste, collapse=".")
  return(overlap)
}

#Function to parallelize
runBoot <- function(antF,beeF,antC,beeC,test){
  while (TRUE){ 
    x <- try(bootOverlap(antF,beeF,antC,beeC,test)) ##Protects from annoying errors that happen with resampling
    if (!inherits(x,"try-error")){
      df <- read.csv(paste(test,"Boot.csv",sep=""))
      df = df[,-c(1)]
      df <- rbind(df,x)
      write.csv(df,file=paste(test,"Boot.csv",sep=""))
    } 
  }
}

runBoot(factorA,factorB,ant,bee,"caste")
runBoot(factorA,factorB,ant,bee,"social")


