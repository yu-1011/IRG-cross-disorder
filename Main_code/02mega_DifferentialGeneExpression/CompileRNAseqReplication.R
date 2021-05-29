#5c_CompileRNAseqReplication

rm(list=ls());options(scipen=999)
library(reshape); library(ggplot2); library(gridExtra)
immunegene <- read.delim("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList.txt",header=T,row.names = 1)

plot_pdf=T

pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  rho = cor(ds1,ds2,method="spearman")
  return(list(slope,intercept, rho))
}

setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")
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

dat= data.frame(SCZ=gvex$SCZ.logFC, ASD=rnaseq.asd$logFC, BD=gvex$BD.logFC, AD=rosmap.ad$AD.logFC, PD=rnaseq.pd$log2FoldChange,dataset="BrainGVEX")
dat= rbind(dat, data.frame(SCZ=cmc$SCZ.logFC, ASD=rnaseq.asd$logFC, BD=cmc$BD.logFC, AD=rosmap.ad$AD.logFC,PD=rnaseq.pd$log2FoldChange, dataset="CommonMind"))
dat2= melt(dat,id=c(1,5))

dat2$value = as.numeric(dat2$value)
colnames(dat2)[3] = "comparison"
fit = pcreg(gvex$SCZ.logFC, rnaseq.asd$logFC); linreg = data.frame(dataset="BrainGVEX", comparison="ASD", slope=fit[[1]], intercept=fit[[2]], color="#F8766D", rho=fit[[3]])
fit = pcreg(gvex$SCZ.logFC, gvex$BD.logFC); linreg = rbind(linreg, data.frame(dataset="BrainGVEX", comparison="BD", slope=fit[[1]], intercept=fit[[2]], color="#00BA38", rho=fit[[3]]))
fit = pcreg(gvex$SCZ.logFC, rosmap.ad$AD.logFC); linreg = rbind(linreg, data.frame(dataset="BrainGVEX", comparison="AD", slope=fit[[1]], intercept=fit[[2]], color="#92683E", rho=fit[[3]]))
fit = pcreg(gvex$SCZ.logFC, rnaseq.pd$log2FoldChange); linreg = rbind(linreg, data.frame(dataset="BrainGVEX", comparison="PD", slope=fit[[1]], intercept=fit[[2]], color="#92683E", rho=fit[[3]]))

fit = pcreg(cmc$SCZ.logFC, rnaseq.asd$logFC); linreg = rbind(linreg,data.frame(dataset="CommonMind", comparison="ASD", slope=fit[[1]], intercept=fit[[2]], color="#F8766D",rho=fit[[3]]))
fit = pcreg(cmc$SCZ.logFC, cmc$BD.logFC); linreg = rbind(linreg,data.frame(dataset="CommonMind", comparison="BD", slope=fit[[1]], intercept=fit[[2]], color="#00BA38",rho=fit[[3]]))
fit = pcreg(cmc$SCZ.logFC, rosmap.ad$AD.logFC); linreg = rbind(linreg,data.frame(dataset="CommonMind", comparison="AD", slope=fit[[1]], intercept=fit[[2]], color="#92683E",rho=fit[[3]]))
fit = pcreg(gvex$SCZ.logFC, rnaseq.pd$log2FoldChange); linreg = rbind(linreg, data.frame(dataset="CommonMind", comparison="PD", slope=fit[[1]], intercept=fit[[2]], color="#92683E", rho=fit[[3]]))

dat2$comparison <- factor(dat2$comparison,levels=c("AD","BD","ASD","PD"))

TxSlope.rnaseq=  ggplot(dat2,aes(x=SCZ,y=value,color=comparison)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(data=linreg, aes(color=comparison, slope=slope, intercept=intercept)) +  facet_grid(.~dataset) +
  scale_colour_manual(values=c("#F8766D","#00BA38","#92683E","#863A76")) +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) + ggtitle("RNAseq Transcriptome Severity") + geom_text(data=linreg, size=3, aes(x=1,y=-1, hjust=1, label=paste("slope=", signif(linreg$slope,2), "\nrho=", signif(linreg$rho,2))),position = position_fill(reverse = T))
if(plot_pdf) pdf(file="../paper_pre/figure/FigS7-RNAseq Replication.pdf", width=10,height=10)
TxSlope.rnaseq
if(plot_pdf) dev.off()



#----Repeat with qSVA normalized data
rnaseq.asd.qsva = read.csv("./results/pretables/RNAseq_ASD_4region_sumstats_qSVA.csv")
gvex.qsva = read.csv("./results/pretables/RNAseq_SCZ_BD_GVEX.csv")
cmc.qsva = read.csv("./results/pretables/RNAseq_SCZ_BD_CMC.csv")

rnaseq.asd.qsva=rnaseq.asd.qsva[match(all_genes,rnaseq.asd.qsva$X),]
gvex.qsva=gvex.qsva[match(all_genes,gvex.qsva$X),]
cmc.qsva = cmc.qsva[match(all_genes, cmc.qsva$X),]

datqsva= data.frame(SCZ=gvex.qsva$SCZ.logFC, ASD=rnaseq.asd.qsva$logFC, BD=gvex.qsva$BD.logFC, dataset="BrainGVEX")
datqsva= rbind(datqsva, data.frame(SCZ=cmc.qsva$SCZ.logFC, ASD=rnaseq.asd.qsva$logFC, BD=cmc.qsva$BD.logFC, dataset="CommonMind"))
datqsva2= melt(datqsva,id=c(1,4))
datqsva2$value = as.numeric(datqsva2$value)
colnames(datqsva2)[3] = "comparison"

fitqsva = pcreg(gvex.qsva$SCZ.logFC, rnaseq.asd.qsva$logFC); linreg.qsva = data.frame(dataset="BrainGVEX", comparison="ASD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#F8766D", rho=fitqsva[[3]])
fitqsva = pcreg(gvex.qsva$SCZ.logFC, gvex.qsva$BD.logFC); linreg.qsva = rbind(linreg.qsva, data.frame(dataset="BrainGVEX", comparison="BD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#00BA38", rho=fitqsva[[3]]))
fitqsva = pcreg(cmc.qsva$SCZ.logFC, rnaseq.asd.qsva$logFC); linreg.qsva = rbind(linreg.qsva,data.frame(dataset="CommonMind", comparison="ASD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#F8766D",rho=fitqsva[[3]]))
fitqsva = pcreg(cmc.qsva$SCZ.logFC, cmc.qsva$BD.logFC); linreg.qsva = rbind(linreg.qsva,data.frame(dataset="CommonMind", comparison="BD", slope=fitqsva[[1]], intercept=fitqsva[[2]], color="#00BA38",rho=fitqsva[[3]]))

TxSlope.rnaseq.qsva=  ggplot(datqsva2,aes(x=SCZ,y=value,color=comparison)) + geom_point(alpha=.5) + 
  geom_abline(slope=1, lty=2) + xlim(-1,1) + ylim(-1,1) + 
  geom_abline(data=linreg.qsva, aes(color=comparison, slope=slope, intercept=intercept)) +  facet_grid(.~dataset) +
  scale_colour_manual(values=c("#F8766D","#00BA38")) +
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) + ggtitle("RNAseq Transcriptome Severity (qSVA)") + 
  geom_text(data=linreg.qsva, size=3, aes(x=1,y=-1, hjust=1, label=paste("slope=", signif(linreg.qsva$slope,2), "\nrho=", signif(linreg.qsva$rho,2))),position = position_fill(reverse = T))






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
datLabel$dataset = factor(datLabel$dataset, levels=c("ROSMAP", "CommonMind", "ASD-pancortical","Dumitriu"))
  
RNAseq.rep=  ggplot(to_plot,aes(x=Microarray,y=RNAseq,color=dataset)) + geom_point(alpha=.5) + 
  geom_smooth(method="lm",fullrange=T) + 
  geom_abline(slope=1, intercept = 0, lty=2) +
  xlab("Microarray log2FC") + ylab("RNAseq log2FC") + coord_fixed(ratio=1) +
  #xlim(c(-1,1)) + ylim(c(-1,1)) +
  ggtitle(" RNAseq Replication")+ facet_grid(.~Dx) +  
  geom_text(data=datLabel, size=3, aes(x=1,y=3,label=paste0("rho=", signif(rho,2))),position = position_fill(reverse = T)) 
            
if(plot_pdf) pdf(file="../paper_pre/figure/FigS7-RNAseq Replication1226.pdf", width=10,height=10)
grid.arrange(grobs=list(TxSlope.rnaseq, RNAseq.rep),ncol=1)
if(plot_pdf) dev.off()



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


