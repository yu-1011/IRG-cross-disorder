#1q_RNAseq_SCZ_ROSMAP.R

rm(list=ls()); options(stringsAsFactors = FALSE)
library(biomaRt); library(WGCNA); library(cqn); library(corrplot); library(ggplot2); library(sva); library(limma)
plot_pdf=T

setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")


datMeta = read.csv("./raw_data/RNAseq_ROSMAP//ROSMAP_clinicalv2.csv")
IDkey = read.csv("./raw_data/RNAseq_ROSMAP/ROSMAP_IDkey.csv")
datMeta = merge(datMeta,IDkey,by = "projid")
datMeta = datMeta[datMeta$mrna_id !="",]
datMeta = datMeta[!duplicated(datMeta$mrna_id),]
rownames(datMeta) = datMeta$mrna_id
datMeta$Batch = as.factor(datMeta$study)
datMeta$Pmi = datMeta$Pmi

Age = datMeta$age_death
Age[Age=="90+"] = 90
Age = as.numeric(Age)
datMeta$Age = Age

datMeta$Dx=factor(datMeta$cogdx)
Dx=factor(datMeta$cogdx)
Dx = as.numeric(Dx)
Dx = ifelse(Dx>=4,"AD","CTL")
datMeta$Dx = factor(Dx, levels=c("CTL", "AD")) 

datExpr = read.table("./raw_data/RNAseq_ROSMAP/ROSMAP_RNAseq_FPKM_gene.tsv",header = F, row.names = 1)
colnames(datExpr) = datExpr[1,]
rownames(datExpr) = gsub("\\..*","",rownames(datExpr))
datExpr = datExpr[-1,]


idx = intersect(rownames(datMeta),colnames(datExpr))
datExpr = datExpr[,idx]
datMeta = datMeta[idx,]
data = apply(datExpr,2,as.numeric)
rownames(data) = rownames(datExpr)
datExpr = data;rm(data)

##------Annotate Probes
getinfo <- c("ensembl_gene_id","external_gene_id","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gc_content")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="feb2014.archive.ensembl.org")
datProbes <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(datExpr),mart=mart)
datProbes = datProbes[match(rownames(datExpr), datProbes$ensembl_gene_id),]
datProbes$length = datProbes$end_position - datProbes$start_position
to_keep = !is.na(datProbes$length)
table(to_keep)
datProbes = datProbes[to_keep,]
datExpr = datExpr[to_keep,]
rownames(datProbes) = datProbes$ensembl_gene_id

to_keep = !is.na(datMeta$Dx)
datMeta = datMeta[to_keep,]; datExpr = datExpr[,to_keep]


#Check consistency
all(colnames(datExpr) %in% rownames(datMeta))

#-----CQN normalize
cqn.dat <- cqn(datExpr,lengths = as.numeric(datProbes$length), x = as.numeric(datProbes$percentage_gc_content),
               lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values
table(is.na(cqn.dat))
datExpr= as.data.frame(cqn.dat)

##-----Filter out genes with low counts: threshold of 1 FPKM in at least 50% of samples
pres = apply(datExpr>1,1,sum) 
to_keep = (pres > 0.5*ncol(datExpr)) ## 15040 genes
table(to_keep)
datExpr = datExpr[to_keep,]
datProbes = datProbes[to_keep,]

save(file="./working_data/RNAseq/RNAseq_AD_ROSMAP.rdata",datExpr,datMeta,datProbes)


if(plot_pdf) pdf("./results/figures/RNAseqQC///SCZ_BD_RNAseq-CommonMind_QC.pdf",width=15,height=8)
par(mfrow=c(3,4), mar=c(2,5,2,2), oma=c(0,0,0,0))

##----------------QC Post-Normalization, Outlier Removal ----------------
## Remove outliers based on network connectivity z-scores
normadj <- (0.5+0.5*bicor(datExpr))^2  ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$Dx,pch=19, main="Outlier Detection", xlab="", ylab="Standardized Network\nConnectivity (Z score)")
legend("bottomright",legend = levels(datMeta$Dx), col = 1:3,pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]; datMeta = datMeta[!outliers,]; 

mod = model.matrix(~datMeta$Dx)
combat = ComBat(datExpr, batch=factor(datMeta$Batch), mod=mod, prior.plots = F)
datExpr = combat


boxplot(datExpr,range=0, col=as.numeric(datMeta$Dx), main="Expression Boxplot)",xaxt = "n")
legend("topright", levels(datMeta$Dx), col=c(1:length(levels(datMeta$Dx))), pch=16, cex=0.8)

plot(density(datExpr[,1]), col=as.numeric(datMeta$Dx)[1], main="Density", ylim=c(0,0.2))
for(i in 2:dim(datExpr)[2]) {
  lines(density(datExpr[,i]), col=as.numeric(datMeta$Dx)[i])  
}
legend("topright", levels(datMeta$Dx), col=c(1:length(levels(datMeta$Dx))), pch=16, cex=0.8)

par(mfrow=c(3,8), mar=c(2,5,2,2), oma=c(0,0,0,0))
#Covariate Plots

    datMeta = datMeta; datExpr = datExpr
    plot_cols = c("black", "red")

  mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  plot(mds$points, col=plot_cols, pch=16, main="MDS: Dx", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));  
  legend("bottomright", levels(datMeta$Dx), col=plot_cols, pch=16, cex=0.8)
  
  plot(datMeta$Dx, ylim=c(0,200),col = plot_cols, main="Subjects")
  A = anova(lm(as.numeric(as.factor(datMeta$msex)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot(as.factor(datMeta$msex) ~ datMeta$Dx, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$Age)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$Age) ~ datMeta$Dx, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$pmi)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$pmi) ~ datMeta$Dx, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$pH)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$pH) ~ datMeta$Dx, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_isolation_RIN)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_isolation_RIN) ~ datMeta$Dx, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
  A = chisq.test(datMeta$Batch,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Batch), col=plot_cols, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
  A = chisq.test(datMeta$Cluster,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Cluster), col=plot_cols, main=paste("Ancestry Cluster,\np=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_Sequencing_Mapped_Reads)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_Sequencing_Mapped_Reads) ~ datMeta$Dx, main=paste("Mapped Reads, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric((datMeta$DLPFC_RNA_Sequencing_Genes_Detected)) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta$DLPFC_RNA_Sequencing_Genes_Detected) ~ datMeta$Dx, main=paste("Genes Detected, p=", signif(p,2)), ylab="", xlab="")
  for(s in paste("seqPC",1:2,sep=""))   {
    A = anova(lm(as.numeric((datMeta[,s])) ~ datMeta$Dx));   p = A$"Pr(>F)"[1];   plot((datMeta[,s]) ~ datMeta$Dx, main=paste(s, ", p=", signif(p,2)), ylab="", xlab="")
  }
  if(length(levels(factor(datMeta$Institution))) > 1) A = chisq.test(datMeta$Institution,datMeta$Dx); p = A$p.value; plot(datMeta$Dx ~ as.factor(datMeta$Institution), col=plot_cols, main=paste("Institution, p=", signif(p,2)), ylab="", xlab="")

if(plot_pdf) dev.off()




#-------Calculate differential expression:  linear model
sumstats = vector(mode="list", length=1); names(sumstats)= c("AD")
mod.ad = model.matrix(~Dx+Age+msex+pmi+Batch,data=datMeta)
datExpr<-datExpr[,rownames(mod.ad)]
datMeta <- datMeta[rownames(mod.ad),]
fit.ad = eBayes(lmFit(datExpr, mod.ad), trend=T,robust=T)

sumstats$AD = topTable(fit.ad, coef=2, number = Inf, sort.by = "none", confint = T)

rnaseq.rosmap = do.call("cbind", sumstats)
write.csv(file="./results/tables/RNAseq_AD_ROSMAP0331.csv", rnaseq.rosmap)
