#3e_network_moduleTrait.R

rm(list=ls()); options(stringsAsFactors = F)
#install.packages("pSI"); library(pSI); library(pSI.data)
library(WGCNA);library(ggplot2); 
library(reshape); library(nlme)

PLOT_PDF = T

rootdir = "/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/"; setwd(rootdir)
load("./working_data//NetworkAnalysis//All_datasets_combined_190410.RData")
load("./working_data//NetworkAnalysis//Brain-Combined-signed-networkConstruction-immunegene0411.RData")


all_colors = unique(moduleColors)
all_colors = all_colors[!grepl("grey",all_colors)]
all_genes = moduleColors
names(all_genes) = rownames(datExpr)

#Step 1 - Calculate module-trait P values, beta, and SEM
datMeta$Group = factor(datMeta$Group,levels = c("CTL","AD", "PD","ETOH","ASD","SCZ","BD","MDD"))
#datMeta$Group = factor(datMeta$Group,levels = c("CTL","AD"))
design=model.matrix(~0+ datMeta$Group)
colnames(design)=levels(datMeta$Group)
moduleColors <- moduleColors
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(datExpr), moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples=dim(datMeta)[1])
moduleTraitCor <- as.data.frame(moduleTraitCor)
moduleTraitCor$AD <- -(moduleTraitCor$AD)
moduleTraitCor <- as.matrix(as.data.frame(moduleTraitCor))

pdf("./results/figures/WGCNA/WGCNA_brainnetworkandtrait0815.pdf",width=10,height=6)

# Will display correlations and their p-values
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
row_idx = which(rownames(moduleTraitCor) %in% c( "MEyellow", "MEcyan","MEtan"))
col_idx = c("ASD","SCZ","BD","MDD","AD","PD")
to_plot = moduleTraitCor[row_idx,col_idx]
textMatrix = paste(signif(moduleTraitCor[row_idx,col_idx], 2), "\n(",
                   signif(moduleTraitPvalue[row_idx,col_idx], 1), ")", sep = "");
dim(textMatrix) = dim(to_plot)
labeledHeatmap(Matrix = to_plot,
               xLabels = colnames(to_plot),
               yLabels = rownames(to_plot),
               ySymbols = rownames(to_plot),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#construct immuneGene network
library(reshape2)
immunegene <- read.table("/Users/normacy/Desktop/immuneGene/01data/ImmuneGeneList",header=T,row.names = 1)
colnames(tom) <- rownames(datExpr)
rownames(tom) <- rownames(datExpr)
net <- melt(tom)
net <- net[net$value>0.1,]
net <- net[net$value!=1,]
write.table(net,file="./working_data/NetworkAnalysis/allnet0816.txt",quote = F,col.names = T,sep = "\t")

idx = which(rownames(datExpr) %in% rownames(immunegene))
itom <- tom[idx,idx]
inet<-melt(itom)
inet<-inet[inet$value>0.1,] #建议对网络的关联值进行一下过滤，过低的就不要了吧

#module是否富集免疫基因
mod <- as.matrix(NA,nrow=length(names(MEs$eigengenes)),ncol=5)
row_idx = which(rownames(moduleTraitCor) %in% c( "MEyellow", "MEcyan","MEtan"))
col_idx = c("ASD","SCZ","BD","MDD","AD","PD")



bpdata = melt(as.matrix(moduleTraitCor[row_idx,col_idx]))
#semdata = melt(moduleTraitSE[row_idx,col_idx])
pdata = melt(moduleTraitPvalue[row_idx,col_idx])
#bpdata$sem = semdata$value
bpdata$p = pdata$value
bpdata$p.symbol = ""
bpdata$p.symbol[bpdata$p<0.1] = "#"
bpdata$p.symbol[bpdata$p<0.05] = "*"
bpdata$p.symbol[bpdata$p<0.01] = "**"
bpdata$p.symbol[bpdata$p<0.001] = "***"

bpdata$X1 = gsub("ME","",bpdata$X1)
bpdata$X1 = factor(bpdata$X1, levels=c("yellow", "tan", "cyan"))

Fig3C.byModule=ggplot(bpdata, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X1,ncol=1, scales="free") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  #geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="beta", x="") + 
#  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  geom_text(color="red",size=4,aes(y=value+ sign(value)*.01), position=position_dodge(.9))  + 
    scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))
Fig3C.byModule
ggsave("./results/figures/WGCNA/WGCNA_brainnetworkandtrait0815v2.pdf",Fig3C.byModule,width = 2.5,height=6)
FigS21.SEMbyModule=ggplot(bpdata, aes(x=X2, y=as.numeric(sem),fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X1,ncol=4, scales="free_x") + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="Std Error of Beta", x="") + 
  #geom_text(color="red",size=4,aes(y=sem*1.1), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))
FigS21.SEMbyModule

Fig3C.byDisease= ggplot(bpdata, aes(x=X1, y=value,fill=X1,group=X1, label=p.symbol))+ facet_wrap(~X2,ncol=1) + 
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(bpdata$X1)) + 
  labs(y="beta", x="") + 
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(angle=45, hjust=1))

if(PLOT_PDF) ggsave(Fig3C.byDisease, filename="./results/figures/Manuscript/Fig3C0809.pdf",width = 2.5,height=10)


dat= bpdata[bpdata$X1 %in% c("yellow","tan", "cyan"),]
neuron=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Neuron Modules") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(neuron, filename="./figures/WGCNA/Modules.Neuron.pdf", width=6,height=3)


dat= bpdata[bpdata$X1 %in% c("yellow"),]
astro=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Astrocyte Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(astro, filename="./figures/WGCNA/Modules.Astro.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("greenyellow"),]
microglia=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Microglial Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(microglia, filename="./figures/WGCNA/Modules.Microglia.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("tan"),]
endo=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Endothelial Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(endo, filename="./figures/WGCNA/Modules.Endo.pdf", width=4,height=3)

dat= bpdata[bpdata$X1 %in% c("blue"),]
bluemod=ggplot(dat, aes(x=X2, y=value,fill=X1,group=X1, label=p.symbol)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +   
  geom_errorbar(aes(ymin=(value - sem), ymax=(value + sem)), position=position_dodge(.9), size=0.3,width=.3) +
  theme_minimal() + scale_fill_manual(name="Group",values=levels(factor(dat$X1))) + 
  labs(y="Differential\nExpression", x="") + ggtitle("Blue Module") +
  geom_text(color="red",size=4,aes(y=value+ sign(value)*sem + sign(value)*.01), position=position_dodge(.9))  + 
  scale_x_discrete() + 
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size=8))
#ggsave(bluemod, filename="./figures/WGCNA/Modules.Blue.pdf", width=4,height=3)

#Module eigengene vs Age trajectory
for (i in 1:dim(MEs)[2]){
  print(names(MEs)[i])
  me = MEs[,i]
  dat=data.frame(Eigengene=me, Age=datMeta$Age, Group=datMeta$Group)
  g5=ggplot(dat,aes(x=Age,y=Eigengene,color=Group)) + geom_point(alpha=0.25) + geom_smooth(span=5, method = "loess")
  ggsave(g5, filename=paste0("./results/figures/WGCNA/Modules.",names(MEs)[i],"Age.pdf"), width=4,height=3)
}
