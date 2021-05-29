##1l) Microarray_IBD_Noble
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

home= "~/Github/" ## insert your GitHub home directory, ie. "C://Users/me/GitHub/"
rootdir = paste(home,"Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap",sep="")
setwd(rootdir)

#1) Load Data
#------------
if(!file.exists("./working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData")) {
  
  datMeta=read.csv("./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_datMeta2.csv")
  rownames(datMeta)= datMeta$GSM
  datMeta$Sex = as.factor(datMeta$Gender)
  datMeta$Ethnicity = as.factor(datMeta$Ethnicity)
  datMeta$Group = as.factor(gsub("Normal", "CTL", gsub("UC", "IBD", datMeta$Disease)))
  datMeta$Tissue = as.factor(datMeta$Anatomic_Location)
  datMeta$Inflammation_State = as.factor(datMeta$Inflammation_State)
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$Study= "Noble.IBD"
  datMeta$Array = "Agilent_G4112A"
  datMeta$Batch = as.factor(datMeta$Run_Date)
  
  ## Read in Expression Data
  filenames=list.files("./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW/")
  RG = read.maimages(files=filenames,source="agilent",path="./raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW")
  plotMD(RG)
  RGb = backgroundCorrect(RG,method= "normexp",offset=50)
  plotDensities(RGb)
  MA <- normalizeWithinArrays(RGb,method="loess")
  plotDensities(MA)
  MA.q <- normalizeBetweenArrays(MA,method="quantile")
  plotDensities(MA.q)
  datExpr = getEAWP(MA.q)$exprs
  rownames(datExpr) = getEAWP(MA.q)$probes$ProbeName
  colnames(datExpr) = gsub(".txt", "", colnames(datExpr))
  
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]
  
  to_remove = which(datMeta$Gender=="unknown")
  datMeta = datMeta[-to_remove,]
  datExpr = datExpr[,-to_remove]
  
  
  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  a = listAttributes(ensembl); f=listFilters(ensembl)
  identifier <- "efg_agilent_wholegenome_4x44k_v1"
  getinfo <- c("efg_agilent_wholegenome_4x44k_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$efg_agilent_wholegenome_4x44k_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  
  save(file="./working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData", datExpr, datMeta, datProbes)
  
  
  #QC Plots before normalization
  datExpr = log2(getEAWP(MA)$exprs)
  par(mfrow=c(2,2))
  boxplot(datExpr, range=0, col = as.numeric(datMeta$Group), main = "Array Boxplot")
  legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
  i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
  for(i in 2:dim(datExpr)[2])
    lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Group[i]))
  legend("topright", legend=levels(datMeta$Group), fill = c(1:length(levels(datMeta$Group))), cex=0.7)
  mds = cmdscale(dist(t(datExpr)))
  plot(mds, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS Plot")
  plot(mds, col=as.numeric(datMeta$Batch), pch=16, main = "MDS by batch")
  
  par(mfrow=c(3,3))
  plot(datMeta$Group, ylim=c(0,150), ylab="Number", main="Subjects")
  A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Sex ~ datMeta$Group, main=paste("Sex, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm((datMeta$Age) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Age ~ datMeta$Group, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Ethnicity) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Ethnicity, main=paste("Ethnicity, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  A = anova(lm(as.numeric(datMeta$Tissue) ~ datMeta$Group)); p = A$"Pr(>F)"[1]
  plot(datMeta$Group ~ datMeta$Tissue, main=paste("Tissue, p=", signif(p,2)), cex.axis = 0.5, ylab="", xlab="")
  
  par(mfrow=c(1,1))
  tree = hclust(dist(t(datExpr)),method="average")
  plot(tree, cex=0.1)
}
