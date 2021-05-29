rm(list=ls()); options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape)
library(NMF); library(WGCNA);library(ggthemes);library(ggrepel)
library(pheatmap);library(corrplot);library(UpSetR);library("clusterProfiler");library("org.Hs.eg.db")
setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")
immunegene <- read.delim("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList.txt",header=T,row.names = 1)
immunegeneFunction <- read.table("/Users/normacy/Desktop/2018immuneGene/01data/AdaptiveImmuneGeneFunction",header=T,sep="\t")

#innate_genes
immunegene <- immunegene[intersect(rownames(immunegene),immunegeneFunction[immunegeneFunction$Catalog=="Adaptive immune response",]$ensembl),]

#brain
asd_meta = read.csv("./results/tables/Microarray_ASD_metaanalysis_190402.csv", row.names=1)
scz_meta = read.csv("./results/tables/Microarray_SCZ_metaanalysis_190402.csv", row.names=1)
bd_meta = read.csv("./results/tables/Microarray_BD_metaanalysis_190402.csv", row.names=1)
mdd_meta = read.csv("./results/tables/Microarray_MDD_metaanalysis_190402.csv", row.names=1)
ad_meta = read.csv("./results/tables/Microarray_ad_metaanalysis_0401.csv", row.names=1)
pd_meta = read.csv("./results/tables/Microarray_PD_metaanalysis_190402.csv", row.names=1)


###
all_genes = intersect(intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(mdd_meta)),rownames(ad_meta)),rownames(pd_meta))
all_genes = intersect(all_genes,rownames(immunegene))
allmeta = matrix(NA,nrow=length(all_genes), 6)
allmeta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
allmeta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
allmeta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
allmeta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
allmeta[,5] = ad_meta$beta[match(all_genes, rownames(ad_meta))]
allmeta[,6] = pd_meta$beta[match(all_genes, rownames(pd_meta))]

#volcano plot
asd_meta <- asd_meta[all_genes,]
ad_meta <- ad_meta[all_genes,]
bd_meta <- bd_meta[all_genes,]
mdd_meta <- mdd_meta[all_genes,]
pd_meta <- pd_meta[all_genes,]
scz_meta <- scz_meta[all_genes,]

asd_meta <- asd_meta[order(asd_meta$fdr),]
ad_meta <- ad_meta[order(ad_meta$fdr),]
bd_meta <- bd_meta[order(bd_meta$fdr),]
mdd_meta <- mdd_meta[order(mdd_meta$fdr),]
pd_meta <- pd_meta[order(pd_meta$fdr),]
scz_meta <- scz_meta[order(scz_meta$fdr),]

asd_meta$disorder <- "ASD"
ad_meta$disorder <- "AD"
bd_meta$disorder <- "BD"
mdd_meta$disorder <- "MDD"
pd_meta$disorder <- "PD"
scz_meta$disorder <- "SCZ"

dfm <- rbind(asd_meta,ad_meta,bd_meta,mdd_meta,pd_meta,scz_meta)
dfm1 <- dfm[dfm$fdr<0.05,]
top5 <- rbind(asd_meta[1:5,],ad_meta[1:5,],bd_meta[1:5,],mdd_meta[1:5,],pd_meta[1:5,],scz_meta[1:5,])

dfm <- na.omit(dfm)
dfm1 <- na.omit(dfm1)
top5 <- na.omit(top5)

ad <- ggplot(dfm[dfm$disorder=="AD",]) +
  geom_point(data = dfm[dfm$disorder=="AD",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="AD",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="AD",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)
asd <- ggplot(dfm[dfm$disorder=="ASD",]) +
  geom_point(data = dfm[dfm$disorder=="ASD",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="ASD",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="ASD",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)
bd <- ggplot(dfm[dfm$disorder=="BD",]) +
  geom_point(data = dfm[dfm$disorder=="BD",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="BD",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="BD",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)
mdd <- ggplot(dfm[dfm$disorder=="MDD",]) +
  geom_point(data = dfm[dfm$disorder=="MDD",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="MDD",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="MDD",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)
pd <- ggplot(dfm[dfm$disorder=="PD",]) +
  geom_point(data = dfm[dfm$disorder=="PD",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="PD",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="PD",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)
scz <- ggplot(dfm[dfm$disorder=="SCZ",]) +
  geom_point(data = dfm[dfm$disorder=="SCZ",],aes(x = beta, y = -log10(fdr)),color = "lightgrey",cex = 0.5) +
  geom_point( data = dfm1[dfm1$disorder=="SCZ",], aes(x = beta, y = -log10(fdr)),color = "darkgrey",cex = 0.5) +
  #geom_point(data = top5, aes(x = beta, y = -log10(fdr)),color = "black", cex = 1.5) +
  theme_bw() +
  xlab("Fold Change") +
  ylab("FDR") +
  #geom_vline(xintercept = 2, col = "red", linetype = "dotted", size = 1) +
  #geom_vline(xintercept = -2,col = "red", linetype = "dotted",size = 1 ) +
  geom_hline(yintercept = -log10(0.05),col = "red", linetype = "dotted", size = 1)+
  geom_label_repel(data = top5[top5$disorder=="SCZ",],aes(x = beta, y = -log10(fdr), label = symbol),cex=3) +
  facet_wrap(~disorder)

valcano <- plot_grid(ad,asd,bd,mdd,pd,scz, align = "h",nrow=2)
valcano + theme(legend.position = 'none',
                axis.text.x = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                panel.grid = element_blank(),
                plot.margin = margin(b = 30))


array.asd = read.csv("./results/tables/Microarray_ASD_metaanalysis_190402.csv")
array.scz = read.csv("./results/tables/Microarray_SCZ_metaanalysis_190402.csv")
array.bd = read.csv("./results/tables/Microarray_BD_metaanalysis_190402.csv")
array.ad = read.csv("./results/tables/Microarray_AD_metaanalysis_0401.csv")
array.pd = read.csv("./results/tables/Microarray_PD_metaanalysis_190402.csv")
array.mdd = read.csv("./results/tables/Microarray_MDD_metaanalysis_190402.csv")

rnaseq.asd = read.csv("./results/pretables/RNAseq_ASD_4region_sumstats_qSVA.csv")
gvex = read.csv("./results/pretables/RNAseq_SCZ_BD_GVEX.csv")
cmc = read.csv("./results/pretables/RNAseq_SCZ_BD_CMC.csv")
rosmap.ad = read.csv("./results/tables/RNAseq_AD0331.csv")
rnaseq.mdd = read.csv("./results/tables/RNAseq_MDD_meta_sumstats.csv")
rnaseq.pd = read.table("./results/tables/RNAseq_PD_sumstats.txt",header=T)


all_genes.array = intersect(array.mdd$X,intersect(array.pd$X,intersect(array.ad$X,intersect(array.asd$X, intersect(array.scz$X, array.bd$X)))))
all_genes.rnaseq = intersect(rnaseq.pd$X,intersect(rosmap.ad$X,intersect(rnaseq.asd$X,intersect(gvex$X,cmc$X))))
all_genes = intersect(all_genes.array, all_genes.rnaseq)
all_genes = intersect(all_genes, rownames(immunegene))

array.asd=array.asd[match(all_genes,array.asd$X),]
array.scz=array.scz[match(all_genes,array.scz$X),]
array.bd=array.bd[match(all_genes,array.bd$X),]
array.ad=array.ad[match(all_genes,array.ad$X),]
array.mdd=array.mdd[match(all_genes,array.mdd$X),]
array.pd=array.pd[match(all_genes,array.pd$X),]

rnaseq.asd=rnaseq.asd[match(all_genes,rnaseq.asd$X),]
gvex=gvex[match(all_genes,gvex$X),]
cmc = cmc[match(all_genes, cmc$X),]
rosmap.ad = rosmap.ad[match(all_genes, rosmap.ad$X),]
rnaseq.pd=rnaseq.pd[match(all_genes,rnaseq.pd$X),]
rnaseq.mdd=rnaseq.mdd[match(all_genes,rnaseq.mdd$X),]

# Individual Datasets
sig.asd = array.asd$fdr<.05
sig.scz = array.scz$fdr<.05
sig.bd = array.bd$fdr<.05
sig.ad= array.ad$fdr<.05
sig.pd= array.pd$fdr<.05
sig.mdd= array.mdd$fdr<.05


to_plot = rbind(data.frame(Microarray=array.scz$beta[sig.scz], RNAseq=gvex$SCZ.logFC[sig.scz], Group="SCZ-GVEX", Dx="SCZ", dataset="BrainGVEX"),
                data.frame(Microarray=array.asd$beta[sig.asd], RNAseq=rnaseq.asd$logFC[sig.asd], Group="ASD", Dx="ASD", dataset='ASD-pancortical'),
                data.frame(Microarray=array.bd$beta[sig.bd], RNAseq=gvex$BD.logFC[sig.bd], Group="BD-GVEX", Dx="BD", dataset="BrainGVEX"),
                data.frame(Microarray=array.bd$beta[sig.bd], RNAseq=cmc$BD.logFC[sig.bd], Group="BD-CMC", Dx="BD", dataset="CommonMind"),
                data.frame(Microarray=array.scz$beta[sig.scz], RNAseq=cmc$SCZ.logFC[sig.scz], Group="SCZ-CMC", Dx="SCZ", dataset="CommonMind"),
                data.frame(Microarray=array.ad$beta[sig.ad], RNAseq=rosmap.ad$AD.logFC[sig.ad], Group="AD-ROSMAP", Dx="AD", dataset="ROSMAP"),
                data.frame(Microarray=array.pd$beta[sig.pd], RNAseq=rnaseq.pd$log2FoldChange[sig.pd], Group="PD-Dumitriu", Dx="PD", dataset="Dumitriu"))
#to_plot$Group = factor(to_plot$Group, levels=c("ASD","SCZ-GVEX", "BD-GVEX", "SCZ-CMC","BD-CMC"))
to_plot$Dx = factor(to_plot$Dx, levels=c("ASD", 'SCZ', 'BD','AD','PD'))

c = cor.test(array.asd$beta[sig.asd], rnaseq.asd$logFC[sig.asd],method="spearman")
datLabel = data.frame(Dx="ASD", dataset='ASD-pancortical', rho=c$estimate, p=c$p.value)
c = cor.test(array.scz$beta[sig.scz], gvex$SCZ.logFC[sig.scz],method="spearman")  #0.85
datLabel = rbind(datLabel,data.frame(Dx="SCZ", dataset='BrainGVEX', rho=c$estimate, p=c$p.value))
c=cor.test(array.bd$beta[sig.bd], gvex$BD.logFC[sig.bd],method="spearman")   #0.78
datLabel = rbind(datLabel,data.frame(Dx="BD", dataset='BrainGVEX', rho=c$estimate, p=c$p.value))
c= cor.test(array.scz$beta[sig.scz], cmc$SCZ.logFC[sig.scz],method="spearman")   #0.65
datLabel = rbind(datLabel,data.frame(Dx="SCZ", dataset='CommonMind', rho=c$estimate, p=c$p.value))
c= cor.test(array.bd$beta[sig.bd], cmc$BD.logFC[sig.bd],method="spearman")   #0.56
datLabel = rbind(datLabel,data.frame(Dx="BD", dataset='CommonMind', rho=c$estimate, p=c$p.value))
c= cor.test(array.ad$beta[sig.ad], rosmap.ad$AD.logFC[sig.ad],method="spearman")   #0.67
datLabel = rbind(datLabel,data.frame(Dx="AD", dataset='ROSMAP', rho=c$estimate, p=c$p.value))
c= cor.test(array.pd$beta[sig.pd], rnaseq.pd$log2FoldChange[sig.pd],method="spearman")   #0.66
datLabel = rbind(datLabel,data.frame(Dx="PD", dataset='Dumitriu', rho=c$estimate, p=c$p.value))

datLabel$Dx = factor(datLabel$Dx, levels=c("ASD", "SCZ", "BD","AD","PD"))
datLabel$dataset = factor(datLabel$dataset, levels=c("ROSMAP","BrainGVEX", "CommonMind", "ASD-pancortical","Dumitriu"))

ad=  ggplot(to_plot[to_plot$Dx=="AD",],aes(y=Microarray,x=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  #geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("RNAseq log2FC") + 
  #coord_fixed(ratio=0.5) +
  scale_colour_manual(values=c("#B274AE")) +
  #xlim(c(-2,2)) + ylim(c(-0.5,0.5)) +
  #ggtitle(" RNAseq Replication")+ 
  facet_grid(~Dx)  + theme(legend.position="none") 
  #geom_text(data=datLabel[datLabel$Dx=="AD",], size=3, aes(x=0,y=2,label=paste0("rho=", signif(rho,2)))) 

scz=  ggplot(to_plot[to_plot$Dx=="SCZ",],aes(x=Microarray,y=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  #geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("") + 
  #coord_fixed(ratio=2) +
  scale_colour_manual(values=c("#9EA021","#29B177")) +
  #xlim(c(-2,2)) + ylim(c(-1,1)) +
  #ggtitle(" RNAseq Replication")+ 
  facet_grid(~Dx)  + theme(legend.position="none") 
  #geom_text(data=datLabel[datLabel$Dx=="SCZ",], size=3, aes(x=1,y=3,label=paste0("rho=", signif(rho,2))),position = position_fill(reverse = T)) 

bd=  ggplot(to_plot[to_plot$Dx=="BD",],aes(x=Microarray,y=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  #geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("") + 
  #coord_fixed(ratio=2) +
  scale_colour_manual(values=c("#9EA021","#29B177")) +
  #xlim(c(-2,2)) + ylim(c(-1,1)) +
  #ggtitle(" RNAseq Replication")+ 
  facet_grid(~Dx)  + theme(legend.position="none") 
  #geom_text(data=datLabel[datLabel$Dx=="BD",], size=3, aes(x=1,y=3,label=paste0("rho=", signif(rho,2))),position = position_fill(reverse = T)) 

asd=  ggplot(to_plot[to_plot$Dx=="ASD",],aes(x=Microarray,y=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  #geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("") + 
  #coord_fixed(ratio=1) +
  scale_colour_manual(values=c("#F8766D","#F8766D")) +
  #xlim(c(-2,2)) + ylim(c(-1,1)) +
  #ggtitle(" RNAseq Replication")+
  facet_grid(~Dx)  + theme(legend.position="none")  
 #geom_text(data=datLabel[datLabel$Dx=="ASD",], size=3, aes(x=1,y=3,label=paste0("rho=", signif(rho,2))),position = position_fill(reverse = T)) 

pd=  ggplot(to_plot[to_plot$Dx=="PD",],aes(y=Microarray,x=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  #geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("") + 
  #coord_fixed(ratio=0.5) +
  scale_colour_manual(values=c("#B274AE")) +
  #xlim(c(-2,2)) + ylim(c(-0.5,0.5)) +
  #ggtitle(" RNAseq Replication")+ 
  facet_grid(~Dx)  + theme(legend.position="none") 
  #geom_text(data=datLabel[datLabel$Dx=="PD",], size=3, aes(x=0.1,y=-0.2,label=paste0("rho=", signif(rho,2)))) 

plot_grid(ad, asd, bd,scz,pd,
          nrow = 1,
          #labels = NULL,
          axis = "l",
          #label_size = 12,
          align = "v"
)

#Inidividual Table
source("./code/00_scripts/fisher_overlap.R")
replicationTable=data.frame()
fisher = data.frame(t(ORA(all_genes[sign(array.asd$beta)==sign(rnaseq.asd$logFC) & rnaseq.asd$P.Value<0.05], all_genes[sig.asd], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="ASD", RNAseq="ASD-pancortical", Array.numDGE=sum(na.omit(sig.asd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.scz$beta)==sign(gvex$SCZ.logFC) & gvex$SCZ.P.Value<0.05], all_genes[sig.scz], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="SCZ", RNAseq="BrainGVEX", Array.numDGE=sum(na.omit(sig.scz)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.bd$beta)==sign(gvex$BD.logFC) & gvex$BD.P.Value<0.05], all_genes[sig.bd], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="BD", RNAseq="BrainGVEX", Array.numDGE=sum(na.omit(sig.bd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))


fisher = data.frame(t(ORA(all_genes[sign(array.scz$beta)==sign(cmc$SCZ.logFC) & cmc$SCZ.P.Value<0.05], all_genes[sig.scz], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="SCZ", RNAseq="CMC", Array.numDGE=sum(na.omit(sig.scz)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.bd$beta)==sign(cmc$BD.logFC) & cmc$BD.P.Value<0.05], all_genes[sig.bd], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="BD", RNAseq="CMC", Array.numDGE=sum(na.omit(sig.bd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.ad$beta)==sign(rosmap.ad$AD.logFC) & rosmap.ad$AD.P.Value<0.05], all_genes[sig.ad], all_genes,all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="AD", RNAseq="ROSMAP", Array.numDGE=sum(na.omit(sig.ad)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.pd$beta)==sign(rnaseq.pd$log2FoldChange) & rnaseq.pd$pvalue<0.05], all_genes[sig.pd], all_genes, all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="PD", RNAseq="Dumitriu", Array.numDGE=sum(na.omit(sig.pd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))

fisher = data.frame(t(ORA(all_genes[sign(array.mdd$beta)==sign(rnaseq.mdd$log2FC) & rnaseq.mdd$pvalue<0.05], all_genes[sig.mdd], all_genes, all_genes)))
replicationTable = rbind(replicationTable, data.frame(Group="MDD", RNAseq="MDD-RNA-seq", Array.numDGE=sum(na.omit(sig.mdd)), 
                                                      RNAseq.concordant=fisher$Overlap, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      PercentOverlap=fisher$X..List.Overlap))


replicationTable

colnames(allmeta) <- c("ASD","SCZ","BD","MDD","AD","PD")
rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)
#allmeta <- na.omit(allmeta)
cor_all<-cor(allmeta,use="pairwise.complete.obs",method="spearman")

##Calculate the slope of transcriptome overlap using principle components regression
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

dat= data.frame(ASD=allmeta$ASD, SCZ=allmeta$SCZ, BD=allmeta$BD, MDD=allmeta$MDD,AD=allmeta$AD,PD=allmeta$PD)
rownames(dat) <- all_genes
dat2= melt(dat,id=2)
dat2$value = as.numeric(dat2$value)
dat2$genes <- rep(all_genes,5)

fit.asd = pcreg(allmeta$SCZ, allmeta$ASD)
fit.bd = pcreg(allmeta$SCZ, allmeta$BD)
fit.mdd = pcreg(allmeta$SCZ, allmeta$MDD)
fit.ad = pcreg(allmeta$SCZ, allmeta$AD)
fit.pd = pcreg(allmeta$SCZ, allmeta$PD)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("^ASD", paste("ASD, r2=", signif(cor_all[2,1],2), sep=""), dat2$variable)
dat2$variable = gsub("^BD", paste("BD, r2=", signif(cor_all[2,3],2), sep=""), dat2$variable)
dat2$variable = gsub("^MDD", paste("MDD, r2=", signif(cor_all[2,4],2), sep=""), dat2$variable)
dat2$variable = gsub("^PD", paste("PD, r2=", signif(cor_all[2,5],2), sep=""), dat2$variable)
dat2$variable = gsub("^AD", paste("AD, r2=", signif(cor_all[2,6],2), sep=""), dat2$variable)

TxSlope.array=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) +
  geom_abline(slope=1, lty=2) + xlim(c(-2,2)) + ylim(c(-2.1,2.1))+
  geom_smooth(method="lm",fullrange=T) + 
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 
pdf("../paper_pre/figure/innate_disorder斜率.pdf")
TxSlope.array
dev.off()

na.allmeta <- na.omit(allmeta)
head(allmeta)
to_plot<-na.allmeta[,c("AD","PD","MDD","ASD","SCZ","BD")]
pdf("../paper_pre/figure/innate_disorderPairsv2.pdf")
chart.Correlation(to_plot)
dev.off()

#cluster
mydata <- na.omit(na.allmeta) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables
# Ward Hierarchical Clustering
d <- dist(t(mydata), method = "minkowski") # distance matrix
fit <- hclust(d, method="single")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

all_beta = matrix(NA,nrow=length(all_genes), 6)
all_beta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
all_beta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
all_beta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
all_beta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
all_beta[,5] = ad_meta$beta[match(all_genes, rownames(ad_meta))]
all_beta[,6] = pd_meta$beta[match(all_genes, rownames(pd_meta))]
all_beta=as.data.frame(all_beta)

all_pvals = matrix(NA,nrow=length(all_genes), 6)
all_pvals[,1] = asd_meta$fdr[match(all_genes, rownames(asd_meta))]
all_pvals[,2] = scz_meta$fdr[match(all_genes, rownames(scz_meta))]
all_pvals[,3] = bd_meta$fdr[match(all_genes, rownames(bd_meta))]
all_pvals[,4] = mdd_meta$fdr[match(all_genes, rownames(mdd_meta))]
all_pvals[,5] = ad_meta$fdr[match(all_genes, rownames(ad_meta))]
all_pvals[,6] = pd_meta$fdr[match(all_genes, rownames(pd_meta))]
all_pvals=  as.data.frame(all_pvals)

colnames(all_beta) = colnames(all_pvals) = c("ASD", "SCZ", "BD", "MDD","AD","PD")
rownames(all_beta) = rownames(all_pvals) =all_genes

all_sig = ifelse(all_pvals<0.05,1,0)
all_sig = as.data.frame(all_sig)
all_sig$genes <-  all_genes
a <- cbind(all_sig,all_beta)
rownames(a) <- all_genes
a <- na.omit(a)
#upset(a, sets = c("ASD", "SCZ", "BD", "MDD", "AAD","AD","PD"),attribute.plots = list(gridrows = 45, plots = list(list(plot = scatter_plot, 
#                                                                                                                     x = "Disease", y = "Gene Number", queries = T), list(plot = scatter_plot, x = "AvgRating", y = "Watches", queries = F)), ncols = 2), query.legend = "bottom")
pdf("../paper_pre/figure/upsetGenes.pdf")
upset(a, sets = c("ASD", "SCZ", "BD", "MDD","AD","PD"), mb.ratio = c(0.55, 0.45))
dev.off()

all_sig$freq <- apply(all_sig[,c(1:6)],1,sum)
atleast_sig <- all_sig[all_sig$freq!=0,]
atleast_beta <- all_beta[atleast_sig$genes,]
atleast_beta <- normalize.quantiles(as.matrix(atleast_beta))
colnames(atleast_beta) = c("ASD", "SCZ", "BD", "MDD","AD","PD")
rownames(atleast_beta) = atleast_sig$genes
pheatmap(t(na.omit(atleast_beta[,c("MDD","ASD","SCZ","BD","AD","PD")])), scale="row",clustering_distance_col = "minkowski", clustering_method = "complete", cluster_row = T)


#Make Corrplot
comparisons = t(combn(seq(1,ncol(allmeta)),2))
barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA)
for (i in 1:dim(comparisons)[1]) {
  x = comparisons[i,1]
  y = comparisons[i,2]
  R = cor.test(allmeta[,x], allmeta[,y])
  rho =cor(allmeta[,x], allmeta[,y], method="spearman", use="pairwise.complete.obs")
  sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allmeta)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allmeta)-3))))/3.92
  barplot[i,] = c(rho, sem, R$p.value)
  rownames(barplot)[i] = paste(colnames(allmeta)[x],colnames(allmeta)[y],sep="-")
}

p_plot <- matrix(NA,ncol=6,nrow=6)
rownames(p_plot) <- c("ASD", "SCZ", "BD", "MDD","AD","PD")
colnames(p_plot) <- rownames(p_plot)
for (m in rownames(p_plot)){
  for (e in colnames(p_plot)) {
    p_plot[m,e] <- cor.test(allmeta[,m],allmeta[,e])$p.value
  }
}
p_plot.fdr = p.adjust(p_plot,method = "fdr")
dim(p_plot.fdr) <- dim(p_plot);rownames(p_plot.fdr) <- c("ASD", "SCZ", "BD", "MDD","AD","PD");colnames(p_plot.fdr) <- rownames(p_plot.fdr)
corrplot(corr <- cor(allmeta,use="pairwise.complete.obs"),order="AOE",type="upper",tl.pos="d", cl.pos="r")
#corrplot(corr = b,add=TRUE, type="lower", method="number",diag=FALSE,tl.pos="n", cl.pos="n")
#corrplot(corr =cor(allmeta,use="pairwise.complete.obs"),p.mat = p_plot.fdr,sig.level = 0.05,insig = "blank",order="AOE",type="upper",tl.pos="tp")
a <- cor(allmeta,use="pairwise.complete.obs")
b <- read.csv("./raw_data/Annotations/genetic-correlation2.csv",header=T,row.names = 1)
col=colorRampPalette(c("navy", "white", "firebrick3"))
corrplot(corr <- a,tl.pos="lt",type = "upper", col=col(10),cl.lim=c(-1,1),addCoef.col="black",p.mat = p_plot.fdr,sig.level = 0.05,insig = "blank",method = "circle",number.cex = 1)
corrplot(corr <- as.matrix(b),type="lower", add = TRUE,cl.pos="n",tl.pos="n" , col=col(10),addCoef.col="black",number.cex =1 )

## Genetic relationship between five psychiatric disorders estimated from genome-wide SNPs. 
## Nat Genet 45, 984–994 (2013).
barplot["SCZ-BD", "genetic.correlation"] = 0.68; barplot["SCZ-BD", "genetic.correlation.SEM"] = 0.04
barplot["SCZ-MDD", "genetic.correlation"] = 0.43; barplot["SCZ-MDD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-SCZ", "genetic.correlation"] = 0.16; barplot["ASD-SCZ", "genetic.correlation.SEM"] = 0.06
barplot["BD-MDD", "genetic.correlation"] = 0.47; barplot["BD-MDD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-BD", "genetic.correlation"] = 0.04; barplot["ASD-BD", "genetic.correlation.SEM"] = 0.06
barplot["ASD-MDD", "genetic.correlation"] = 0.05; barplot["ASD-MDD", "genetic.correlation.SEM"] = 0.09
barplot["ASD-AD", "genetic.correlation"] = 0.0494; barplot["ASD-AD", "genetic.correlation.SEM"] = 0.09
barplot["BD-AD", "genetic.correlation"] = -.02; barplot["BD-AD", "genetic.correlation.SEM"] = 0.05
barplot["AD-PD", "genetic.correlation"] = -.07; barplot["AD-PD", "genetic.correlation.SEM"] = 0.09
barplot["MDD-AD", "genetic.correlation"] = 0.04; barplot["MDD-AD", "genetic.correlation.SEM"] = 0.06
barplot["SCZ-AD", "genetic.correlation"] = 0.03; barplot["SCZ-AD", "genetic.correlation.SEM"] = 0.04
barplot["SCZ-PD", "genetic.correlation"] = 0.05; barplot["SCZ-PD", "genetic.correlation.SEM"] = 0.05
barplot["ASD-PD", "genetic.correlation"] = -.19; barplot["ASD-PD", "genetic.correlation.SEM"] = 0.1
barplot["BD-PD", "genetic.correlation"] = 0.07; barplot["BD-PD", "genetic.correlation.SEM"] = 0.05
barplot["MDD-PD", "genetic.correlation"] = 0.07; barplot["MDD-PD", "genetic.correlation.SEM"] = 0.05

barplot$p.fdr = p.adjust(barplot$p.fdr,method="fdr")
barplot$p.bootstrap = null$prob[findInterval(barplot$Mean, null$null)]
barplot$p.symbol = ""
barplot$p.symbol[barplot$p.fdr<0.05] = "*"
barplot$p.symbol[barplot$p.fdr<0.01] = "**"
barplot$p.symbol[barplot$p.fdr<0.001] = "***"
barplot$Comparison = rownames(barplot)
barplot$modality="innate-genes"

c = cor.test(barplot$Mean, barplot$genetic.correlation, method="pearson", use="pairwise.complete.obs")
CI.95 = CIrho(c$estimate, N=13)[2:3]
e = expression(paste(rho, "=0.46, p.value=0.08"))
linearMod =  lm(formula = barplot$Mean ~ barplot$genetic.correlation)


plot2 = ggplot(barplot , aes(x=genetic.correlation,y=Mean,label=Comparison, color=modality)) + 
  geom_point(size=3) +
  labs(x="Genetic (SNP-based) correlation", y=expression(paste("Transcriptome correlation (", rho, ")", sep=""))) +       
  geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM),colour=modality), width=0,size=0.8) +
  geom_errorbarh(aes(xmin=(genetic.correlation - genetic.correlation.SEM), xmax=(genetic.correlation + genetic.correlation.SEM)), height=0,size=0.8) + 
  theme(
    axis.text.x=element_text(size=10, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,16,2),"mm")
  ) + 
  #geom_text(size=4,aes(colour=genes),position = position_jitter(width = 0, height = .05)) + 
  stat_smooth(method="lm",alpha=0, fullrange=T) +
  geom_label_repel(aes(x=barplot$genetic.correlation,y=barplot$Mean,label=Comparison),
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+
  geom_abline(intercept=linearMod$coefficients[1], slope=linearMod$coefficients[2],lty=2) +
  theme(legend.position="none")+
  annotate("text", label=as.character(e), x=Inf,y=Inf, vjust=1,hjust=1, parse=T)

pdf("../paper_pre/figure/innate_transcriptome&geneticCorr.pdf")
plot2
dev.off()


