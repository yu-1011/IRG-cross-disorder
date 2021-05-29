rm(list=ls()); options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R")
library(ggplot2); library(mada); library(reshape2)
library(NMF); library(WGCNA);library(ggthemes);library(ggrepel)
library(pheatmap);library(corrplot);library(UpSetR);library("clusterProfiler");library("org.Hs.eg.db")
setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")
immunegene <- read.delim("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList.txt",header=T,row.names = 1)
immunegeneFunction <- read.table("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneFunction",header=T,sep="\t")


#brain
asd_meta = read.csv("./results/tables/Microarray_ASD_metaanalysis_190402.csv", row.names=1)
scz_meta = read.csv("./results/tables/Microarray_SCZ_metaanalysis_190402.csv", row.names=1)
bd_meta = read.csv("./results/tables/Microarray_BD_metaanalysis_190402.csv", row.names=1)
mdd_meta = read.csv("./results/tables/Microarray_MDD_metaanalysis_190402.csv", row.names=1)
ad_meta = read.csv("./results/tables/Microarray_ad_metaanalysis_0401.csv", row.names=1)
pd_meta = read.csv("./results/tables/Microarray_PD_metaanalysis_190402.csv", row.names=1)


###
all_genes = intersect(intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(mdd_meta)),rownames(ad_meta)),rownames(pd_meta))
all_genes = intersect(all_genes,rownames(na.omit(mdd_meta)))
all_genes = intersect(all_genes,rownames(immunegene))
allmeta = matrix(NA,nrow=length(all_genes), 6)
allmeta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
allmeta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
allmeta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
allmeta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
allmeta[,5] = ad_meta$beta[match(all_genes, rownames(ad_meta))]
allmeta[,6] = pd_meta$beta[match(all_genes, rownames(pd_meta))]

colnames(allmeta) <- c("ASD","SCZ","BD","MDD","AD","PD")

rownames(allmeta) = all_genes
allmeta=  as.data.frame(allmeta)
#allmeta <- na.omit(allmeta)
cor_all<-cor(allmeta,use="pairwise.complete.obs",method="spearman")
tree= hclust(as.dist(1-cor_all),method="single")
plot(tree)
##Calculate the slope of transcriptome overlap using principle components regression
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  return(list(slope,intercept))
}

dat= data.frame(ASD=allmeta$ASD, SCZ=allmeta$SCZ, BD=allmeta$BD, MDD=allmeta$MDD,AD=allmeta$AD,PD=allmeta$PD)
dat<- as.data.frame(dat)
colnames(dat) <- c("ASD","SCZ","BD","MDD","AD","PD")
dat2= melt(dat,id=2)
dat2$value = as.numeric(dat2$value)
dat2$genes <- rep(all_genes,5)


fit.asd = pcreg(allmeta$SCZ, allmeta$ASD)
fit.bd = pcreg(allmeta$SCZ, allmeta$BD)
fit.mdd = pcreg(allmeta$SCZ, allmeta$MDD)
fit.ad = pcreg(allmeta$SCZ, allmeta$AD)
fit.pd = pcreg(allmeta$SCZ, allmeta$PD)

dat2$variable = as.character(dat2$variable)
dat2$variable = gsub("^ASD", paste("ASD, slope=", signif(fit.asd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("^BD", paste("BD, slope=", signif(fit.bd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("^MDD", paste("MDD, slope=", signif(fit.mdd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("^PD", paste("PD, slope=", signif(fit.pd[[1]],2), sep=""), dat2$variable)
dat2$variable = gsub("^AD", paste("AD, slope=", signif(fit.ad[[1]],2), sep=""), dat2$variable)

TxSlope.array=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) +
  geom_abline(slope=1, lty=2) + xlim(c(-2,2)) + ylim(c(-2.1,2.1))+ geom_point(alpha = 0.1)+
  geom_smooth(method="lm",fullrange=T) + 
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 
pdf("../paper_pre/figure/disorder斜率0419.pdf")
TxSlope.array
dev.off()

na.allmeta <- na.omit(allmeta)
head(allmeta)
library(PerformanceAnalytics)
pdf("../paper_pre/figure/disorderPairsv2.pdf")
to_plot<-na.allmeta[,c("MDD","ASD","SCZ","BD","AD","PD")]
chart.Correlation(to_plot)
dev.off()

#TAC1
TxSlope.TAC1=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) +
  geom_abline(slope=1, lty=2) + xlim(c(-2,2)) + ylim(c(-2.1,2.1))+
  geom_smooth(method="lm",fullrange=T) + 
  geom_point(data = dat2[dat2$genes=="ENSG00000006128",],mapping=aes(SCZ,value),show.legend = T) +
  geom_label_repel(data=dat2[dat2$genes=="ENSG00000006128",],aes(x=dat2[dat2$genes=="ENSG00000006128",]$SCZ,y=dat2[dat2$genes=="ENSG00000006128",]$value,label=c("TAC1")),
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 
pdf("../paper_pre/figure/disorderTAC1.pdf")
TxSlope.TAC1
dev.off()

#TAC1
TxSlope.TAC1=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) +
  geom_abline(slope=1, lty=2) + xlim(c(-2,2)) + ylim(c(-2.1,2.1))+
  geom_smooth(method="lm",fullrange=T) + 
  geom_point(data = dat2[dat2$genes=="ENSG00000006128",],mapping=aes(SCZ,value),show.legend = T) +
  geom_label_repel(data=dat2[dat2$genes=="ENSG00000006128",],aes(x=dat2[dat2$genes=="ENSG00000006128",]$SCZ,y=dat2[dat2$genes=="ENSG00000006128",]$value,label=c("TAC1")),
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 
pdf("../paper_pre/figure/disorderTAC1.pdf")
TxSlope.TAC1
dev.off()


TAC1 <- melt(dat["ENSG00000006128",])
TAC1$p.symbol=c("**","***","*","","***","*")
TAC1$gene<-"TAC1"

TAC1.change = ggplot(TAC1,aes(x = variable, y=value,fill=variable,label=p.symbol))  +  
  geom_bar(stat="identity",width=0.75) +
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm"))+ xlab("Disorders")+ylab("Expression Differences(beta)")+
  geom_text(color="red",size=4,aes(y=value))
pdf("../paper_pre/figure/disorderTAC1change.pdf",width=3,height=6)
TAC1.change
dev.off()

#IL17
IL17 <- melt(dat["ENSG00000112115",])
IL17$p.symbol=c("**","","","*","","")

IL17.change = ggplot(IL17,aes(x = variable, y=value,fill=variable,label=p.symbol))  +  
  geom_bar(stat="identity",width=0.75) +
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm"))+ xlab("Disorders")+ylab("Effect Size")+
  geom_text(color="red",size=4,aes(y=value))
pdf("../paper_pre/figure/disorderIL17change.pdf",width=5,height=5)
IL17.change
dev.off()



TxSlope.CRH=ggplot(dat2,aes(x=SCZ,y=value,color=variable)) +
  geom_abline(slope=1, lty=2) + xlim(c(-2,2)) + ylim(c(-2.1,2.1))+
  geom_smooth(method="lm",fullrange=T) + 
  geom_point(data = dat2[dat2$genes=="ENSG00000147571",],mapping=aes(SCZ,value),show.legend = T) +
  geom_label_repel(data=dat2[dat2$genes=="ENSG00000147571",],aes(x=dat2[dat2$genes=="ENSG00000147571",]$SCZ,y=dat2[dat2$genes=="ENSG00000147571",]$value,label=c("CRH")),
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+
  xlab("Schizophrenia (log2FC)") + ylab("Disease2 (log2FC)") +
  coord_fixed(ratio=1) 
pdf("../paper_pre/figure/disorderCRH.pdf")
TxSlope.CRH
dev.off()

CRH <- melt(dat["ENSG00000147571",])
CRH$p.symbol=c("***","***","**","*","**","")
CRH$gene <- "CRH"

CRH.change = ggplot(CRH,aes(x = variable, y=value,fill=variable,label=p.symbol))  +  
  geom_bar(stat="identity",width=0.75) +
    theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm"))+ xlab("Disorders")+ylab("Expression Differences(beta)")+
    geom_text(color="red",size=4,aes(y=value+value*0.1))
pdf("../paper_pre/figure/disorderCRHchange.pdf")
CRH.change
dev.off()
gene <- rbind(TAC1,CRH)
gene$variable <- factor(gene$variable,levels = c("AD","ASD","BD","MDD","PD","SCZ"))
gene.change = ggplot(gene,aes(x = variable, y=value,fill=variable,label=p.symbol))  +  
  geom_bar(stat="identity",width=0.75) +
  theme(
    axis.text.x=element_text(angle=50, size=12, hjust=1),
    axis.text.y=element_text(size=10, vjust=0.5), 
    legend.title = element_text(size=12), 
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
    axis.title.y = element_text(size=12, vjust=0.5),
    plot.margin=unit(c(2,2,1,2),"mm"))+ xlab("Disorders")+ylab("Effect Size")+
  geom_text(color="red",size=4,aes(y=value+value*0.1))+facet_grid(~gene)

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
corrplot(corr <- a,tl.pos="lt",type = "upper", col=col(10),cl.lim=c(-1,1),addCoef.col="black",p.mat = p_plot.fdr,sig.level = 0.05,insig = "blank",method = "circle",number.cex = 0.0001)
corrplot(corr <- as.matrix(b),type="lower", add = TRUE,cl.pos="n",tl.pos="n" , col=col(10),addCoef.col="black",number.cex =0.0000001 )

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
barplot$modality="immune-genes"

c = cor.test(barplot$Mean, barplot$genetic.correlation, method="pearson", use="pairwise.complete.obs")
CI.95 = CIrho(c$estimate, N=13)[2:3]
e = expression(paste(rho, "=0.47, p.value=0.07"))
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

psychiatric.barplot <- barplot[c("ASD-SCZ","ASD-BD","ASD-MDD","ASD-AD","SCZ-BD","SCZ-MDD","SCZ-AD","BD-MDD","BD-AD","MDD-AD"),]
c = cor.test(psychiatric.barplot$Mean, psychiatric.barplot$genetic.correlation, method="spearman")
CI.95 = CIrho(c$estimate, N=13)[2:3]
e = expression(paste(rho, "=0.66, p.value=0.03"))
linearMod =  lm(formula = psychiatric.barplot$Mean ~ psychiatric.barplot$genetic.correlation)


plot2 = ggplot(psychiatric.barplot , aes(x=genetic.correlation,y=Mean,label=Comparison, color=modality)) + 
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
  geom_label_repel(aes(x=psychiatric.barplot$genetic.correlation,y=psychiatric.barplot$Mean,label=Comparison),
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+
  geom_abline(intercept=linearMod$coefficients[1], slope=linearMod$coefficients[2],lty=2) +
  theme(legend.position="none")+
  annotate("text", label=as.character(e), x=Inf,y=Inf, vjust=1,hjust=1, parse=T)

pdf("../paper_pre/figure/transcriptome&geneticCorr0401.pdf")
plot2
dev.off()

library(ggplot2)
library(gProfileR)
library(gridExtra)
library(grid)

#aad_immune = subset(aad_meta[match(all_genes, rownames(aad_meta)),],fdr<0.05)
ad_immune = subset(ad_meta[match(all_genes, rownames(ad_meta)),],fdr<0.05)
asd_immune = subset(asd_meta[match(all_genes, rownames(asd_meta)),],fdr<0.05)
bd_immune = subset(bd_meta[match(all_genes, rownames(bd_meta)),],fdr<0.05)
#ibd_immune = subset(ibd_meta[match(all_genes, rownames(ibd_meta)),],fdr<0.05)
mdd_immune = subset(mdd_meta[match(all_genes, rownames(mdd_meta)),],fdr<0.05)
pd_immune = subset(pd_meta[match(all_genes, rownames(pd_meta)),],fdr<0.05)
scz_immune = subset(scz_meta[match(all_genes, rownames(scz_meta)),],fdr<0.05) 


#Biological theme comparison
gcpare = list()
gcpare$asd = rownames(subset(asd_meta[match(all_genes, rownames(asd_meta)),],fdr<0.05))
gcpare$scz = rownames(subset(scz_meta[match(all_genes, rownames(scz_meta)),],fdr<0.05))
gcpare$bd = rownames(subset(bd_meta[match(all_genes, rownames(bd_meta)),],fdr<0.05))
gcpare$mdd = rownames(subset(mdd_meta[match(all_genes, rownames(mdd_meta)),],fdr<0.05))
#gcpare$aad = rownames(subset(aad_meta[match(all_genes, rownames(aad_meta)),],fdr<0.05))
#gcpare$ibd= rownames(subset(ibd_meta[match(all_genes, rownames(ibd_meta)),],fdr<0.05))
gcpare$ad = rownames(subset(ad_meta[match(all_genes, rownames(ad_meta)),],fdr<0.05))
gcpare$pd = rownames(subset(pd_meta[match(all_genes, rownames(pd_meta)),],fdr<0.05))

ck <- compareCluster(geneCluster = gcpare, fun = "enrichGO",
                     OrgDb         = org.Hs.eg.db,
                     keytype       = 'ENSEMBL',
                     ont           = 'BP',
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
sim.ck.0.8 <- clusterProfiler::simplify(ck, cutoff=0.8,by = "p.adjust",
                                   measure = "Wang", semData = NULL)

write.table(as.data.frame(ck@compareClusterResult),file="./results/tables/DifferentialGeneAnnotation_200126.txt",sep="\t",row.names = F,col.names = T,quote = F)
formagma <- ck@compareClusterResult[,c("Cluster","ID","geneID")]
formagma <- na.omit(melt(separate(formagma,geneID,sep="/",into=paste("gene",c(1:106))),c("Cluster","ID"),paste("gene",c(1:106))))
adformagma <- formagma[formagma$Cluster=="ad",]
asdformagma <- formagma[formagma$Cluster=="asd",]
bdformagma <- formagma[formagma$Cluster=="bd",]
mddformagma <- formagma[formagma$Cluster=="mdd",]
sczformagma <- formagma[formagma$Cluster=="scz",]
pdformagma <- formagma[formagma$Cluster=="pd",]
write.table(adformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/adpathformagma",quote=F,row.names = F,col.names = T)
write.table(bdformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/bdpathformagma",quote=F,row.names = F,col.names = T)
write.table(mddformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/mddpathformagma",quote=F,row.names = F,col.names = T)
write.table(asdformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/asdpathformagma",quote=F,row.names = F,col.names = T)
write.table(pdformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/pdpathformagma",quote=F,row.names = F,col.names = T)
write.table(sczformagma[,c(2,4)],"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/sczpathformagma",quote=F,row.names = F,col.names = T)

write.table(ag,"/Users/normacy/Desktop/2018immuneGene/paper_pre/table/IRGsformagma",quote=F,row.names = F,col.names = T)



all_ck <- ck@compareClusterResult[which(ck@compareClusterResult$ID%in%c("GO:0046651","GO:0031349","GO:0040017","GO:0051047","GO:2001233","GO:0045088")),]

library(reshape2)
test <-  dcast(sim.ck@compareClusterResult,Description~Cluster,value.var = "qvalue")
rownames(test) <- test[,1]
test <- test[,-1]
test[!is.na(test)] <- 1
test[is.na(test)] <- 0
colnames(test) <- c("ASD","SCZ","BD","MDD","AD","PD")
pdf("../paper_pre/figure/upsetPathway1206.pdf")
upset(test,sets = c("ASD","SCZ","BD","MDD","AD","PD"))
dev.off()
write.csv(test,file="../paper_pre/table/pathway1206.csv")
test$sum <- apply(test,1,sum)
test[test$sum=="6",]

ggplot(all_ck,aes(x=Cluster,y=-log10(qvalue),size=Count,col=Cluster))+geom_point()+facet_wrap(~Description,ncol=3)

data <- sim.ck_bp@compareClusterResult
data <- data[data$Description=="positive regulation of locomotion",]
data$Description <-  factor(data$Description,levels=data$Description)
ggplot(data = data, mapping = aes(x = Cluster , y = -log10(qvalue),fill=Cluster)) +
  ylab("-log10(q.value)")+ geom_bar(stat = 'identity') + 
  geom_hline(yintercept=-log10(0.05),col="red")+
  xlab("Positive Regulation of Locomotion")+
#  coord_flip()+
  theme(axis.title.x = element_text(size = 15),axis.text.y.left = element_text(size=15))
#geom_text(aes(label = Description, vjust = 1.1, hjust = -0.5, angle = 45))+
#facet_wrap(~Type)
#ggplot2::scale_y_discrete(labels=function(x) str_wrap(x, width=55))+ggplot2::scale_size(range=c(2, 6))


write.table(as.data.frame(ck_bp),file="./results/tables/DifferentialGeneAnnotation_bp_191120.txt",sep="\t",row.names = F,col.names = T,quote = F)


mydf <- data.frame(ID=rownames(ad_immune),FC=ad_immune$beta)
mydf$group<- "NA"
mydf$group <- ifelse(mydf$FC>0,yes="upregulated",no="downregulated")
mydf$othergroup <- "ad"

mydf1 <- data.frame(ID=rownames(asd_immune),FC=asd_immune$beta)
mydf1$group<- "NA"
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "asd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(bd_immune),FC=bd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "bd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(mdd_immune),FC=mdd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "mdd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(pd_immune),FC=pd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "pd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(scz_immune),FC=scz_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "scz"
mydf <- rbind(mydf,mydf1)

formula_res_bp <- compareCluster(ID~group+othergroup, data=mydf, fun="enrichGO",
                              OrgDb         = org.Hs.eg.db,
                              keytype       = 'ENSEMBL',
                              ont           = 'BP',
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)


gene.symbols= asd_meta$symbol[match(all_genes, rownames(asd_meta))]

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
library( preprocessCore )
atleast_beta <- normalize.quantiles(as.matrix(atleast_beta))
colnames(atleast_beta) = c("ASD", "SCZ", "BD", "MDD","AD","PD")
rownames(atleast_beta) = atleast_sig$genes

fc <- immunegeneFunction[immunegeneFunction$ensembl %in% rownames(atleast_beta),]
annotation_col = matrix(NA, nrow=length(rownames(atleast_beta)), ncol=length(unique(fc$Type)))
rownames(annotation_col) = rownames(atleast_beta)
colnames(annotation_col) = unique(fc$Type)

for (i in colnames(annotation_col)){
  for (m in c(1:length(rownames(annotation_col)))){
  annotation_col[m,i] = ifelse(rownames(annotation_col)[m] %in% fc$ensembl[fc$Type==i],i,"na")
  }}

annotation <- annotation_col[rownames(annotation_col) %in% colnames(t(na.omit(atleast_beta))),]
annotation <- annotation[rownames(annotation) %in% rownames(atleast_beta),]
atleast_beta <- atleast_beta[rownames(atleast_beta) %in% rownames(annotation),]
rownames(annotation) <- 1:nrow(annotation)
rownames(atleast_beta) <- 1:nrow(atleast_beta)
pheatmap(t(atleast_beta))

pheatmap(normalizeQuantiles(allmeta))
quantno(allmeta[,1])

write.table(a,file="../paper_pre/table/a.txt",quote=F,col.names = T)

data <- read.table(pipe("pbpaste"),sep='\t',header=T)
data$fraction = data$n / sum(data$n)
data$ymax = cumsum(data$fraction)
data$ymin = c(0, head(data$ymax, n = -1))
data$lb= paste(data$Type,"\n(",signif(data$fraction,2)*100,"% )")

p <- ggplot(data = data, aes(fill = Type, ymax = ymax, ymin = ymin, xmax = 4, xmin = 3)) +
  geom_rect(show.legend = F,alpha=0.8) +
  scale_fill_brewer(palette = 'Set3')+
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "") + 
  xlim(c(0, 5)) +
  theme_light() +
  theme(panel.grid=element_blank()) + ## 去掉白色外框
  theme(axis.text=element_blank()) + ## 把图旁边的标签去掉
  theme(axis.ticks=element_blank()) + ## 去掉左上角的坐标刻度线
  theme(panel.border=element_blank()) + ## 去掉最外层的正方形边框
  #geom_label_repel(aes(x = 4, y = ((ymin+ymax)/2),label = Description) ,size=4)+
  geom_label_repel(aes(x = 2, y = (ymax+ymin)/2,label = lb) ,size=4)
 
pdf("../paper_pre/figure/DiffPathway.pdf")
p
dev.off()
p


mydf <- data.frame(ID=rownames(ad_immune),FC=ad_immune$beta)
mydf$group<- "NA"
mydf$group <- ifelse(mydf$FC>0,yes="upregulated",no="downregulated")
mydf$othergroup <- "ad"

mydf1 <- data.frame(ID=rownames(asd_immune),FC=asd_immune$beta)
mydf1$group<- "NA"
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "asd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(bd_immune),FC=bd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "bd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(mdd_immune),FC=mdd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "mdd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(pd_immune),FC=pd_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "pd"
mydf <- rbind(mydf,mydf1)

mydf1 <- data.frame(ID=rownames(scz_immune),FC=scz_immune$beta)
mydf1$group <- ifelse(mydf1$FC>0,yes="upregulated",no="downregulated")
mydf1$othergroup <- "scz"
mydf <- rbind(mydf,mydf1)

formula_res <- compareCluster(ID~group+othergroup, data=mydf, fun="enrichGO",
                              OrgDb         = org.Hs.eg.db,
                              keytype       = 'ENSEMBL',
                              ont           = 'BP',
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

