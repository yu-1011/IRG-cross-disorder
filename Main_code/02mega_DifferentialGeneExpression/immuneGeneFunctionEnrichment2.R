#3g_Annotation_Pathway.R

rm(list=ls()); options(stringsAsFactors = F)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA)


#load("./working_data//NetworkAnalysis//brainCombinedfinalizedNetwork_0809.RData")
immunegene <- read.table("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList",header=T,sep="\t")
immunegeneFunction <- read.table("/Users/normacy/Desktop/immuneGene/01data/ImmuneGeneFunction",header=T,sep="\t")

#
asd_meta = read.csv("./results/tables//Microarray_ASD_immuneGENE_metaanalysis_062718.csv", row.names=1)
scz_meta = read.csv("./results/tables/Microarray_SCZ_immunegene_metaanalysis_062718.csv", row.names=1)
bd_meta = read.csv("./results/tables/Microarray_BD_immunegene_metaanalysis_062718.csv", row.names=1)
mdd_meta = read.csv("./results/tables/Microarray_MDD_immunegene_metaanalysis_062718.csv", row.names=1)
aad_meta = read.csv("./results/tables/Microarray_AAD_immunegene_metaanalysis_062718.csv", row.names=1)
ibd_meta = read.csv("./results/tables//Microarray_IBD_immunegene_metaanalysis_062718.csv", row.names=1)
ad_meta = read.csv("./results/tables//Microarray_AD_immunegene_metaanalysis_062718.csv", row.names=1)
pd_meta = read.csv("./results/tables//Microarray_PD_immunegene_metaanalysis_092017.csv", row.names=1)

asd_meta$threshold <- as.factor(ifelse (asd_meta$fdr < 0.05  ,ifelse( asd_meta$beta >0 ,'Up','Down'),'Not'))

aad_meta$threshold <- as.factor(ifelse (aad_meta$fdr < 0.05  ,ifelse( aad_meta$beta >0 ,'Up','Down'),'Not'))
ad_meta$threshold <- as.factor(ifelse (ad_meta$fdr < 0.05  ,ifelse( ad_meta$beta >0 ,'Up','Down'),'Not'))
bd_meta$threshold <- as.factor(ifelse (bd_meta$fdr < 0.05  ,ifelse( bd_meta$beta >0 ,'Up','Down'),'Not'))
ibd_meta$threshold <- as.factor(ifelse (ibd_meta$fdr < 0.05 ,ifelse( ibd_meta$beta >0 ,'Up','Down'),'Not'))
mdd_meta$threshold <- as.factor(ifelse (mdd_meta$fdr < 0.05 ,ifelse( mdd_meta$beta >0 ,'Up','Down'),'Not'))
pd_meta$threshold <- as.factor(ifelse (pd_meta$fdr < 0.05 ,ifelse( pd_meta$beta >0 ,'Up','Down'),'Not'))
scz_meta$threshold <- as.factor(ifelse (scz_meta$fdr < 0.05 ,ifelse( scz_meta$beta >0 ,'Up','Down'),'Not'))

  
  
source("./code/00_scripts/fisher_overlap.R")
table.p = matrix(NA, nrow=8, ncol=5)
num.p = table.p
rownames(table.p) = c("asd","ad","aad","bd","ibd","mdd","pd","scz"); rownames(num.p) = c("asd","aad","ad","bd","ibd","mdd","pd","scz")
colnames(table.p) = paste("M", 1:5,sep="");colnames(num.p) = paste("M", 1:5,sep="")
table.or = table.p; num.p = table.p; table.ov = table.p
#hgnc = multiExpr[[1]]$datProbes$`Ensembl Gene ID`

#rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Up"),])

#q = length(intersect(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Up"),]),immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e]))
#b = length(intersect(rownames(get(paste0(m,"_meta"))), immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e]))
#n = length(rownames(get(paste0(m,"_meta")))) - b
#k = length(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Up"),]))
#phyper(q,b,n,k)

for (m in rownames(table.p)) {
  for (e in colnames(table.p)){
    f = ORA(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Up"),]),immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e], rownames(get(paste0(m,"_meta"))), immunegeneFunction$ensembl)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
    table.ov[m,e] = as.numeric(f[[5]])
  }
}

for (m in rownames(table.p)) {
  for (e in colnames(table.p)){
    q = length(intersect(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Down"),]),immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e]))
    b = length(intersect(rownames(get(paste0(m,"_meta"))), immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e]))
    n = length(rownames(get(paste0(m,"_meta")))) - b
    k = length(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Down"),]))
    table.p[m,e] = phyper(q,b,n,k)
    table.ov[m,e] = q
  }
}

for (m in rownames(table.p)) {
  for (e in colnames(table.p)){
    f = ORA(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Down"),]),immunegeneFunction$ensembl[immunegeneFunction$Module.Label==e], rownames(get(paste0(m,"_meta"))), immunegeneFunction$ensembl)
    table.or[m,e] = as.numeric(f[[1]])
    table.p[m,e] = as.numeric(f[[2]])
    table.ov[m,e] = as.numeric(f[[5]])
  }
}

for (m in rownames(table.p)) {
  for (e in colnames(table.p)){
    num.p[m,e] = paste(signif(table.p[m,e],1),"\n(",table.ov[m,e],")")
  }
}

table.p.fdr = p.adjust(table.p,method ="fdr")
dim(table.p.fdr) = dim(table.p); dimnames(table.p.fdr) = dimnames(table.p)


to_plot = table.p.fdr
head(to_plot)
row <- {}
col <- {}
a <- {}
for (m in rownames(to_plot)) {
   a <- (paste(m,"\n(",length(rownames(get(paste0(m,"_meta"))[which(get(paste0(m,"_meta"))$threshold=="Up"),])),")"))
   row <- append(row,a)
}
rownames(to_plot) <- row

for (m in colnames(to_plot)) {
  a <- (paste(m,"(",length(immunegeneFunction$ensembl[immunegeneFunction$Module.Label==m]),")"))
  col <- append(col,a)
}
colnames(to_plot)<-col


heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000], 
          scale="none",trace="none",cexRow = 1.0,cexCol = 1.0, density.info = "none",
          colsep=0:19,rowsep=0:13,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5,  keysize=1.1, xlab="Immune Gene Function",
          key=T,key.xlab="-log10(P)", cellnote=num.p, notecex=1.2, notecol="black",main="Immune Gene Function\nNetwork Enrichment")

dev.off()



