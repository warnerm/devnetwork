load("~/GitHub/devnetwork/data/DEtests.RData")
library(cowplot)
DE_pallete = c("#ffffcc","#a1dab4","#41b6c4","#225ea8")

main_theme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1, color = "black",fill = NA),
        text=element_text(family='Arial'),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 22))

############
###Fig 1a
############
#Compare DE definition at each stage
antRes[[2]]$ortholog_found = antRes[[2]]$OGG_found = FALSE
antRes[[2]]$ortholog_found[antRes[[2]]$Gene %in% AllPS$Gene.x] = TRUE
antRes[[2]]$OGG_found[antRes[[2]]$Gene %in% ACUogg$gene_Mphar] = TRUE
aM = melt(antRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
beeRes[[2]]$ortholog_found = beeRes[[2]]$OGG_found = FALSE
beeRes[[2]]$ortholog_found[beeRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeRes[[2]]$OGG_found[beeRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

aD$DEboth = bD$DEboth = allD$DEboth
aD$OGG = aD$OGG - aD$DEboth
bD$OGG = bD$OGG - allD$DEboth

aDM = melt(aD,id.vars = "variable")
bDM = melt(bD,id.vars = "variable")
colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")
aDM$species = "ant"
bDM$species = "bee"

d = rbind(aDM,bDM)
d$species=as.factor(d$species)
levels(d$species) = c("M. pharaonis","A. mellifera")
levels(d$DEtype) = c("NDE","no ortholog","paralogs present","ortholog present,\nnon-conserved caste-bias","ortholog present,\nconserved caste-bias")

p1 <- ggplot(d[d$DEtype!="NDE",],
       aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  ylim(0,4500)+
  facet_grid(. ~ species)+
  scale_fill_manual(values = DE_pallete)+
  ylab("number of caste-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="lightgrey"),
        plot.margin = margin(2,2,2,2,"cm"))

p2 <- ggplot(d[d$DEtype!="NDE" & d$species=="M. pharaonis",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  scale_fill_manual(values = DE_pallete)+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  ylab("proportion of DEGs")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))


p3 <- ggplot(d[d$DEtype!="NDE" & d$species=="A. mellifera",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  scale_fill_manual(values = DE_pallete)+
  ylab("proportion of DEGs")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))

png("~/GitHub/devnetwork/figures/Fig1a.png",height=2500,width=3500,res=300)
ggdraw()+
  draw_plot(p1+
              theme(legend.text = element_text(size=13),
                    legend.key.width = unit(1,"cm")))+
  draw_plot(p2,x=0.2,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.58,height=0.18,width=0.18)
dev.off()

############
###Fig 1b
############
#Compare DE definition at each stage
antSocRes[[2]]$ortholog_found = antSocRes[[2]]$OGG_found = FALSE
antSocRes[[2]]$ortholog_found[antSocRes[[2]]$Gene %in% AllPS$Gene.x] = TRUE
antSocRes[[2]]$OGG_found[antSocRes[[2]]$Gene %in% ACUogg$gene_Mphar] = TRUE
aM = melt(antSocRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
aD = ddply(aM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
beeSocRes[[2]]$ortholog_found = beeSocRes[[2]]$OGG_found = FALSE
beeSocRes[[2]]$ortholog_found[beeSocRes[[2]]$Gene %in% AllPS$Gene.y] = TRUE
beeSocRes[[2]]$OGG_found[beeSocRes[[2]]$Gene %in% ACUogg$gene_Amel] = TRUE
bM = melt(beeSocRes[[2]],id.vars = c("Gene","ortholog_found","OGG_found"))
bD = ddply(bM,~variable,summarize,
           NDE = sum(value=="nonDE"),
           no_ortholog = sum(value!="nonDE" & !ortholog_found),
           dup = sum(value!="nonDE" & ortholog_found & !OGG_found),
           OGG = sum(value!="nonDE" & OGG_found))
colnames(bM)[5] = "value_apis"
aM = merge(aM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Mphar")
bM = merge(bM[,-c(2,3)], ACUogg,by.x="Gene",by.y="gene_Amel")
allM = merge(aM,bM,by=c("OGGacu","variable"))

allD = ddply(allM,~variable,summarize,
             DEboth = sum(value_apis!="nonDE" & value != "nonDE"))

aD$DEboth = bD$DEboth = allD$DEboth
aD$OGG = aD$OGG - aD$DEboth
bD$OGG = bD$OGG - allD$DEboth

aDM = melt(aD,id.vars = "variable")
bDM = melt(bD,id.vars = "variable")
colnames(aDM) = colnames(bDM) = c("stage","DEtype","value")
aDM$species = "ant"
bDM$species = "bee"

d = rbind(aDM,bDM)
d$species=as.factor(d$species)
levels(d$species) = c("M. pharaonis","A. mellifera")
levels(d$DEtype) = c("NDE","no ortholog","paralogs present","ortholog present,\nnon-conserved role-bias","ortholog present,\nconserved role-bias")
levels(d$stage) = c("head","thorax","abdomen")

p1 <- ggplot(d[d$DEtype!="NDE",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity")+
  main_theme+
  xlab("tissue")+
  ylim(0,4500)+
  facet_grid(. ~ species)+
  scale_fill_manual(values = DE_pallete)+
  ylab("number of role-biased genes")+
  theme(axis.text.x = element_text(angle=-25,hjust=0.1),
        legend.position = "top",
        strip.text = element_text(size=20,face="bold.italic"),
        axis.title.y=element_text(margin=margin(t=0,l=0,r=10,b=0)),
        legend.title = element_blank(),
        strip.background = element_rect(color="black",fill="lightgrey"),
        plot.margin = margin(2,2,2,2,"cm"))

p2 <- ggplot(d[d$DEtype!="NDE" & d$species=="M. pharaonis",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  scale_fill_manual(values = DE_pallete)+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  ylab("proportion of DEGs")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))


p3 <- ggplot(d[d$DEtype!="NDE" & d$species=="A. mellifera",],
             aes(x = stage, y = value, fill = DEtype))+
  geom_bar(stat="identity",position = "fill")+
  main_theme+
  theme(axis.text = element_text(size = 8.5),axis.ticks = element_blank())+
  xlab("")+
  scale_fill_manual(values = DE_pallete)+
  ylab("proportion of DEGs")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_rect(fill="transparent"),
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,r=0,l=0,b=0)),
        axis.text.y = element_text(margin=margin(t=0,r=-3,l=0,b=0)),
        axis.title = element_text(size=10),
        plot.margin = margin(0,0,0,0,"cm"))

png("~/GitHub/devnetwork/figures/Fig1b.png",height=2500,width=3500,res=300)
ggdraw()+
  draw_plot(p1+
              theme(legend.text = element_text(size=13),
                    legend.key.width = unit(1,"cm")))+
  draw_plot(p2,x=0.2,y=0.58,height=0.18,width=0.18)+
  draw_plot(p3,x=0.59,y=0.58,height=0.18,width=0.18)
dev.off()
