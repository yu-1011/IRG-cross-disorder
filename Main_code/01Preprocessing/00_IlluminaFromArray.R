##1A) Microarray_ASD_voineagu

rm(list=ls()); options(stringsAsFactors=F)
suppressPackageStartupMessages(T)
#source("http://bioconductor.org/biocLite.R")


library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(gplots); library(Cairo); library(GEOquery)
library(biomaRt); library(sva)

home= "~/Github/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)


# Step 1) Download and normalize raw microarray data
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_ASD_voineagu_normalized.RData")) {
  #-----------Load Raw Data----------------------------------------
  if(!file.exists("./raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28521&format=file&file=GSE28521%5Fnon%2Dnormalized%5Fdata%2Etxt%2Egz",
                destfile="./raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")}

  data.lumi = lumiR("./raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")
  datMeta = read.csv("./raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_metaData.csv")
  matchSN = match(sampleNames(data.lumi), datMeta$GEO_SampleName)
  datMeta = datMeta[matchSN,]
  
  # log2 transform and create separate datasets for: all samples, cortex, cerebellum
  dataAll.lumi<-lumiT(data.lumi, method="log2"); 
  dataCTX.lumi<-dataAll.lumi[,datMeta$Brain.area!="C"]; 
  dataCBL.lumi<-dataAll.lumi[,datMeta$Brain.area=="C"]
  
  #Normalize
  dataAll_N.lumi<-lumiN(dataAll.lumi, method="quantile"); 
  dataCTX_N.lumi<- lumiN(dataCTX.lumi, method="quantile"); 
  dataCBL_N.lumi<- lumiN(dataCBL.lumi, method="quantile");
  
  #Extract expression data for Cortex Only
  datExpr = exprs(dataCTX_N.lumi)
  datMeta = datMeta[datMeta$Brain.area!="C",]
  datExpr.prenorm = exprs(dataCTX.lumi)
  
  ##Re-annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  identifier <- "illumina_humanref_8_v3"
  getinfo <- c("illumina_humanref_8_v3", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat[,"illumina_humanref_8_v3"])
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  rownames(datProbes) = datProbes[,1]
  
  #Clean Meta-Data
  datMeta$PMI[is.na(datMeta$PMI)]=mean(datMeta$PMI, na.rm=T)
  datMeta$A.C = factor(datMeta$A.C, levels = c("C", "A"))
  datMeta$Brain.area = factor(datMeta$Brain.area, levels = c("F", "T"))
  datMeta$Chip = factor(datMeta$Chip)
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group[datMeta$A.C=="A"] = "ASD"
  datMeta$Group[datMeta$A.C=="C"] = "CTL"
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "ASD"))
  datMeta$Study = "ASD.voineagu"
  
  #Raw data QC Plots
  par(mfrow=c(2,2))
  boxplot(datExpr.prenorm, col = as.numeric(datMeta$A.C),range=0, main="Expression pre-Normalization")
  
  i = 1; plot(density((datExpr.prenorm[,i]), na.rm=T), col = as.numeric(datMeta$A.C)[i], main="Hist of Log2 Exp", xlab = "log2 exp", ylim=c(0,0.8))
  for(i in 2:dim(datExpr.prenorm)[2])
    lines(density((datExpr.prenorm[,i]), na.rm=T), col = as.numeric(datMeta$A.C)[i])
  
  mds = cmdscale(dist(t(datExpr.prenorm)))
  plot(mds, col=as.numeric(as.factor(datMeta$A.C)), pch=16, main="MDS Plot")
  legend("topright", c("CTL", "ASD"), col=c(1:2), pch=16, cex=0.8)
  
  plot(mds, col=as.numeric(as.factor(datMeta$Brain.area)), pch=16, main="MDS Plot")
  legend("topright", c("Frontal", "Temporal"), col=c(1:2), pch=16, cex=0.8)
  
  #Assess initial levels of counfounding
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,30), ylab="Number", main="Samples")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  datMeta$Chip = factor(datMeta$Chip)
  A = anova(lm(as.numeric(datMeta$Chip) ~ datMeta$A.C)); p = A$"Pr(>F)"[1]
  plot(datMeta$A.C ~ datMeta$Chip, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Brain.area) ~ datMeta$A.C)); p = A$"Pr(>F)"[1]
  plot(datMeta$Brain.area ~ datMeta$A.C, main=paste("Brain Region, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$A.C)); p = A$"Pr(>F)"[1]
  plot(datMeta$RIN ~ datMeta$A.C, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$A.C)); p = A$"Pr(>F)"[1]
  plot(datMeta$PMI ~ datMeta$A.C, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_ASD_voineagu_normalized.RData", datExpr, datMeta, datProbes)
  
} 
