#Merge phylostrata results
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
Aps$psName[Aps$psName %in% insect_arthropod] = "insect"
Aps$psName[Aps$psName %in% aculeata] = "aculeata"
Aps$psName[Aps$psName=="Monomorium"] = "novel"
Aps$psName = factor(Aps$psName,levels = c("old","insect","Hymenoptera","aculeata","ant","novel"))
levels(Aps$psName)[3] = "hymenoptera"

Btax = "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Holometabola; Hymenoptera; Apocrita; Aculeata; Apoidea; Apidae; Apinae; Apini; Apis"
Btax = strsplit(Btax,"; ")[[1]]
Bps <- mergeRes2(7460,"Amel",Btax,"novel","Amel_7460_blastAll_ps")
Bps$psName[Bps$psName %in% old] = "old"
Bps$psName[Bps$psName %in% insect_arthropod] = "insect"
Bps$psName[Bps$psName %in% aculeata] = "aculeata"
Bps$psName[Bps$psName %in% bee_name] = "bee"
Bps$psName[Bps$psName == "Apis"] = "novel"
Bps$psName = factor(Bps$psName,levels = c("old","insect","Hymenoptera","aculeata","bee","novel"))
levels(Bps$psName)[3] = "hymenoptera"

AllPS = merge(Aps[!is.na(Aps$OGGacu),],Bps[!is.na(Bps$OGGacu),],by = "OGGacu")
AllPS1 = AllPS[,c(1,4,8)]
AllPS_sum = ddply(AllPS1,~OGGacu,summarize,
                  ps = min(ps.x,ps.y))
AllPS_sum$psName = "old"
AllPS_sum$psName[AllPS$ps > 8] = "insect"
AllPS_sum$psName[AllPS$ps > 18] = "hymenoptera"
AllPS_sum$psName[AllPS$ps > 20] = "aculeata"
AllPS_sum$psName = factor(AllPS_sum$psName,levels = c("old","insect","hymenoptera","aculeata"))

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

#Sex-bias
design <- model.matrix(~sex+time_factor,data=droplevels(f_sex[f_sex$time>25,]))
d = counts[,colnames(counts) %in% f_sex$SRA[f_sex$time>25]]
sexGenes <- EdgeR(d,design,2)
sexGenes$Gene = rownames(sexGenes)
sexGenes = merge(sexGenes,ENDogg,by.x="Gene",by.y="gene_Dmel")

#adding Jasper data
#Calculate tau (tissue specificity) for honey bee genes
calculatetau <- function(factor,expr){
  meanExpr <- lapply(levels(factor$tissue),function(x) 
    rowSums(as.data.frame(expr[,colnames(expr) %in% factor$SRA[factor$tissue==x]]))/sum(factor$tissue==x))
  meanExpr <- as.data.frame(do.call(cbind,meanExpr))
  colnames(meanExpr) = levels(factor$tissue)
  geneDeviation = apply(meanExpr,1,function(x){
    sum(ldply(lapply(c(1:12),function(i) 1-x[i]/max(x))))
  })
  
  tau = geneDeviation/(11)
  
  return(list(tau,meanExpr))
}

counts <- read.csv("~/GitHub/devnetwork/data/Jasper_counts.csv")
fpkm <- read.csv("~/GitHub/devnetwork/data/Jasper_fpkm.csv")
samples <- read.table("~/GitHub/devnetwork/data/samples_jasper.txt",head=T)
rownames(counts) = counts$gene_id
counts = counts[,-c(1)]
rownames(fpkm) = fpkm$gene_id
fpkm = fpkm[,-c(1)]

samples = samples[samples$caste!="callow",]
counts = counts[,colnames(counts) %in% samples$SRA]
beeTissue <- calculatetau(samples,fpkm)
tau <- data.frame(Gene=names(beeTissue[[1]]),tau=beeTissue[[1]])

save(ACUogg,ENDogg,Aps,Bps,AllPS_sum,sexGenes,tau,file="~/GitHub/devnetwork/phylo_results/collectedPhylo.RData")

