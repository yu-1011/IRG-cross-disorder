##1D) Microarray_SCZ_BD_chen
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= '~/Github/' ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

set.seed(100)

#1) Load Data from .cel
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData")) {
  
  datMeta=read.csv("./raw_data/Microarray/Chen_GSE35978/GSE35978_datMeta.csv")
  rownames(datMeta)= datMeta$Chip
  datMeta$Region = as.factor(datMeta$Region)
  datMeta$Group[datMeta$Group=="Bipolar"] = "BD"; datMeta$Group[datMeta$Group=="Control"] = "CTL";  datMeta$Group[datMeta$Group=="Schizophrenia"] = "SCZ";   datMeta$Group[datMeta$Group=="Depression"] = "MDD";
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "MDD", "BD", "SCZ"))
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Sex=as.factor(datMeta$Sex)
  datMeta$Study="SCZ.BD.Chen"

  #Remove NA group
  to_keep = !is.na(datMeta$Group)
  datMeta = datMeta[to_keep,]
  
  ## Get parietal cortex samples only
  ctx_only = rownames(datMeta)[datMeta$Region=="PCTX"]
  all_samples = list.files("./raw_data/Microarray/Chen_GSE35978/GSE35978_CEL/")
  ctx_sample_files = all_samples[pmatch(ctx_only, all_samples)]
  
  ## Read in Expression Data
  data.affy = ReadAffy(celfile.path="./raw_data/Microarray/Chen_GSE35978/GSE35978_CEL", filenames=ctx_sample_files)
  datExpr = affy::rma(data.affy,verbose=T,normalize=T,background=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substring(colnames(datExpr), 1, 9)
  
  ## Get RNA degradation
  RNAdeg = AffyRNAdeg(data.affy)
  plotAffyRNAdeg(RNAdeg)
  RNAdeg$sample.names = substring(RNAdeg$sample.names,1,9)
  idx=  match(rownames(datMeta), RNAdeg$sample.names)
  datMeta$RNAdeg = RNAdeg$slope[idx]
  
  ## Get batch information
  batch = as.factor(substring(protocolData(data.affy)$ScanDate,1,10))
  idx = match(rownames(datMeta), colnames(datExpr))
  datMeta$Batch = batch[idx]
  
  ## Align expression and phenoData matricies
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  a = listAttributes(ensembl)
  identifier <- "affy_hugene_1_0_st_v1"
  getinfo <- c("affy_hugene_1_0_st_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$affy_hugene_1_0_st_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC Plots before normalization
  datExpr = log2(exprs(data.affy))
  par(mfrow=c(2,2))
  boxplot(datExpr, range=0, col = as.numeric(datMeta$Group), main = "Array Boxplot")
  legend("topright", legend=c("CTL", "Depression", "Bipolar", "Schizophrenia"), fill = c(1:4), cex=0.7)
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
  legend("topright", legend=c("CTL", "Depression", "Bipolar", "Schizophrenia"), fill = c(1:4), cex=0.7)
  mds = cmdscale(dist(t(datExpr)))
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")
  
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,50), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$PMI ~ datMeta$Group, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$pH) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$pH ~ datMeta$Group, main=paste("pH, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RNAdeg) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$RNAdeg ~ datMeta$Group, main=paste("RNAdeg, p=", signif(p,2)), ylab="", xlab="")
  
  par(mfrow=c(1,1))
  tree = hclust(dist(t(datExpr)),method="average")
  plot(tree, cex=0.1)
}

