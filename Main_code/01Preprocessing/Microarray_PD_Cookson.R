##1b_Microarray_PD_Cookson.R

rm(list=ls()); options(stringsAsFactors = F)

library(limma); library(WGCNA); library(biomaRt); library(sva)
library(ggfortify)
setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")

#Load & clean metaData
datMeta = read.csv("./raw_data/Microarray/PD_Cookson/datMeta_GSE28894.csv",header = T)
rownames(datMeta) = datMeta$Source
datMeta$disease = factor(datMeta$disease,levels =  c("control", "Parkinson's_disease"))
datMeta$gender = as.factor(datMeta$gender)
datMeta$age = as.numeric(datMeta$age)
datMeta$Sample_source_name = as.factor(datMeta$Sample_source_name)
datMeta$Term_Accession_Number = as.factor(datMeta$Term_Accession_Number)
datMeta$Study="PD.Cookson"

#Load quantile normaliezed micorarray data
datExpr = read.delim(file="./raw_data/Microarray/PD_Cookson/GSE28894_series_matrix.txt",comment="!")
rownames(datExpr) = datExpr[,1]
datExpr=datExpr[,-1]
datExpr <- log(datExpr,base=2)
datExpr <- na.omit(datExpr)
datExpr <- normalizeQuantiles(datExpr)


idx = match(colnames(datExpr), rownames(datMeta))
datMeta = datMeta[idx,]

#Initial QC Plots
boxplot(datExpr, range = 0, col= as.numeric(datMeta$disease), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");  legend("top", legend=levels(datMeta$disease),pch=19,col=c(1:2))
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$disease)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$disease), col=c(1:length(levels(datMeta$disease))), pch=16, cex=0.8)
plot(hclust(dist(t(datExpr)),method="average"))


##Annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
identifier <- "illumina_humanref_8_v3"
getinfo <- c("illumina_humanref_8_v3", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
geneDat<- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
geneDat<-getBM(attributes=c("illumina_humanref_8_v3", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position"), filters="illumina_humanref_8_v3", values=rownames(datExpr),
               mart = ensembl)
#geneDat <- read.delim("./raw_data/Microarray/PD_Cookson/GPL6104-11576.txt",comment.char = "#",header = T)

idx = match(rownames(datExpr), geneDat$ID)
datProbes = cbind(rownames(datExpr), geneDat[idx,])

colnames(datExpr) = rownames(datMeta)
save(file="./working_data/Microarray/01_Normalized/Microarray_PD_Cookson_normalized.RData", datMeta, datExpr, datProbes)
pdf("./results/figures/MicroarrayQC/PD_Cookson_QC.pdf", width=11,height=8.5)
par(mfrow=c(3,4))
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$disease), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
legend("bottomleft",pch=16, legend = levels(datMeta$disease), col = 1:2)
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

model = model.matrix(~disease+Sample_source_name+Term_Accession_Number+age+gender, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_PD_Cookson_normalized_balanced.RData", datMeta, datExpr, datProbes,model)

#Post-Normalization QC
boxplot(datExpr, range = 0, col= as.numeric(datMeta$disease), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   
legend("topright", legend=levels(datMeta$disease), col = 1:2, pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease[i]),) } ;   legend("topright", levels(datMeta$disease), cex=0.7, text.col = 1:2)

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$disease)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$disease), col=unique(as.numeric(datMeta$disease)), pch=16, cex=0.8)

# Plot potential Colinearities
plot(datMeta$disease, ylab="Number", main="Subjects")
#A = anova(lm(as.numeric(datMeta$Term_Accession_Number.1) ~ datMeta$disease)); p = A$"Pr(>F)"[1];   plot(datMeta$disease ~ datMeta$datMeta$Term_Accession_Number.1, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$disease) ~ as.factor(datMeta$Term_Accession_Number))); p = A$"Pr(>F)"[1];   plot(datMeta$disease ~ as.factor(datMeta$Term_Accession_Number), main=paste("Batch.1, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$disease) ~ as.factor(datMeta$Sample_source_name))); p = A$"Pr(>F)"[1];   plot(datMeta$disease ~ as.factor(datMeta$Sample_source_name), main=paste("Brain Regions, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$gender) ~ datMeta$disease)); p = A$"Pr(>F)"[1];   plot(datMeta$gender ~ datMeta$disease, main=paste("gender, p=1"), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$age) ~ datMeta$disease)); p = A$"Pr(>F)"[1];   plot(datMeta$age ~ datMeta$disease, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$Batch) ~ datMeta$disease)); p = A$"Pr(>F)"[1];   plot(datMeta$Batch ~ datMeta$disease, main=paste("Batch.2, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$disease)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$disease, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");

# Cluster Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
gender_col = as.numeric((datMeta$gender)); gender_col[gender_col==2] = "blue"; gender_col[gender_col==1]="pink"
age_col = numbers2colors(datMeta$age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$age),max(datMeta$age)))
disease_col = as.numeric((datMeta$disease)); disease_col[disease_col==2] = "blue";disease_col[disease_col==1]="pink"
plotDendroAndColors(tree, colors = cbind(disease_col, as.numeric(datMeta$Term_Accession_Number.1), gender_col, age_col), groupLabels = c("disease", "gender", "Age"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()

#PCA for sample
pca <-prcomp(t(datExpr),scale=F)
group <-factor(datMeta$disease)
colour_group <-cm.colors(length(unique(group)))
colour <-colour_group[as.numeric(factor(group))]
colour
group2 <-data.frame(group)
data2 <-cbind(t(datExpr),group2)
autoplot(pca,data=data2,colour='group')

# CollapseRows Probes --> Genes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
geneAnno<-getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position"), filters="entrezgene", values=datProbes$Entrez_Gene_ID,
                mart = ensembl)
names(datProbes) <- gsub("Entrez_Gene_ID","entrezgene",names(datProbes))
datProbes <- merge(geneAnno,datProbes,by="entrezgene")
datProbes <- datProbes[!duplicated(datProbes$Search_Key),]
rownames(geneDat) <- geneDat$ensembl_gene_id
realGenes = !is.na(geneDat$ensembl_gene_id)  
table(realGenes)
datExpr1 = datExpr[realGenes,]; datProbes1 = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = geneDat$ensembl_gene_id,rowID = rownames(geneDat))
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], rownames(datProbes))
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Batch Correction
plot(datMeta$disease ~ datMeta$Term_Accession_Number, main = "Batch Balance", ylab ="", xlab = "Batch")
mod = model.matrix(~disease, data=datMeta)
batch = factor(datMeta$Term_Accession_Number)
datExpr = ComBat(datExpr, batch=batch, mod=mod)
dev.off()

#QC Plots
boxplot(datExpr, range=2, col = as.numeric(datMeta$disease), main = "Array Boxplot")
legend("topright", legend=levels(datMeta$disease), fill = c(1:length(levels(datMeta$disease))), cex=0.7)

# Histogram
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease[i]), main="Hist of Log2 Exp", xlab = "log2 exp")
for(i in 2:dim(datExpr)[2])
  lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease[i]))
legend("topright", legend=levels(datMeta$disease), fill = c(1:length(levels(datMeta$disease))), cex=0.7)


## Regress all Technical and Non-DX biological covariates
#datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
X = model.matrix(~disease+Sample_source_name+age+gender, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,3:4]) %*% (as.matrix(beta[3:4,]))) 
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_PD_Cookson_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)
