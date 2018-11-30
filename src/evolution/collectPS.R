#Merge phylostrata results
mergeRes2 <- function(tax_id,species,taxNames,fullName,psFile){
  TGmap <- read.table(paste("../../results/phylostratigraphy/TGmap_",species,".txt",sep=""))
  blast <- read.table(paste("../../results/phylostratigraphy/",tax_id,"_",species,"_blastRes.tab",sep=""))
  PS <- read.table(paste("../../results/phylostratigraphy/",psFile,sep=""))
  ACUogg <- read.table(paste("../../results/phylostratigraphy/",species,"_",tax_id,"_ACU",sep=""))
  
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
  oggGs$ps[oggGs$ps==40]=maxPS #Dummy number indicating novel gene
  oggGs$psName = apply(oggGs,1,function(x) taxNames[as.integer(as.character(x[4]))])
  
  return(oggGs)
}
#merge endopterygota results
mergeEndo <- function(species,num){
  Dmel_END = read.table(paste("../../results/phylostratigraphy/",species,"_",num,"_END",sep=""))
  Dmel_blast =  read.table(paste("../../results/phylostratigraphy/",num,"_",species,"_blastRes.tab",sep=""))
  Dmel_tg = read.table(paste("../../results/phylostratigraphy/TGmap_",species,".txt",sep=""))
  
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
aculeata = c("Aculeata")
hymenoptera = c("Hymenoptera","Apocrita")
bee_name = c("Apoidea", "Apidae", "Apinae","Apini")
Aps$psName[Aps$psName %in% ant_name] = "ant"
Aps$psName[Aps$psName %in% old] = "old"
Aps$psName[Aps$psName %in% insect_arthropod] = "insect"
Aps$psName[Aps$psName %in% aculeata] = "aculeata"
Aps$psName[Aps$psName=="Monomorium"] = "novel"
Aps$psName[Aps$psName %in% hymenoptera] = "hymenoptera"
Aps$psName = factor(Aps$psName,levels = c("old","insect","hymenoptera","aculeata","ant","novel"))
Btax = "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Holometabola; Hymenoptera; Apocrita; Aculeata; Apoidea; Apidae; Apinae; Apini; Apis"
Btax = strsplit(Btax,"; ")[[1]]
Bps <- mergeRes2(7460,"Amel",Btax,"novel","Amel_7460_blastAll_ps")
Bps$psName[Bps$psName %in% old] = "old"
Bps$psName[Bps$psName %in% insect_arthropod] = "insect"
Bps$psName[Bps$psName %in% aculeata] = "aculeata"
Bps$psName[Bps$psName %in% bee_name] = "bee"
Bps$psName[Bps$psName == "Apis"] = "novel"
Bps$psName[Bps$psName %in% hymenoptera] = "hymenoptera"
Bps$psName = factor(Bps$psName,levels = c("old","insect","hymenoptera","aculeata","bee","novel"))
#Merge phylostrata results and get orthologs
AllPS = merge(Aps[!is.na(Aps$OGGacu),],Bps[!is.na(Bps$OGGacu),],by = "OGGacu")
#Take the minumum ps (i.e. oldest)
AllPS_sum = ddply(AllPS,~OGGacu,summarize,
                  ps = min(ps.x,ps.y))
AllPS_sum$psName = "old"
AllPS_sum$psName[AllPS_sum$ps > 8] = "insect"
AllPS_sum$psName[AllPS_sum$ps > 18] = "hymenoptera"
AllPS_sum$psName[AllPS_sum$ps > 20] = "aculeata"
AllPS_sum$psName = factor(AllPS_sum$psName,levels = c("old","insect","hymenoptera","aculeata"))
#Get endopterygota orthologs
AmelEND = mergeEndo("Amel",7460)
MpharEND = mergeEndo("Mphar",307658)
DmelEND = mergeEndo("Dmel",7227)
ENDogg = merge(AmelEND,MpharEND,by="OGGend")
ENDogg = merge(ENDogg,DmelEND,by="OGGend")
colnames(ENDogg) = c("OGGend","gene_Amel","gene_Mphar","gene_Dmel")
#Get 1-1 orthologs for apis and monomorium
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
save(ACUogg,ENDogg,Aps,Bps,AllPS,AllPS_sum,file="../../results/collectedPhylo.RData")