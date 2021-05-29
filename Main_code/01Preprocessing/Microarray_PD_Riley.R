##1b_Microarray_PD_Riley.R

rm(list=ls()); options(stringsAsFactors = F)

library(limma); library(WGCNA); library(biomaRt); library(sva)
library(ggfortify)
setwd("/Users/normacy/Desktop/immuneGene/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master")

#Load & clean metaData
datMeta = read.table("./raw_data/Microarray/PD_Riley/datMeta_GSE54282.csv",header = T)
rownames(datMeta) = datMeta$Assay_Name
datMeta$disease_state = factor(datMeta$disease_state,levels =  c("normal", "PD"))
datMeta$Sex = as.factor(datMeta$sex)
datMeta$age = as.numeric(datMeta$age)
datMeta$Sample_source_name = as.factor(datMeta$Sample_source_name)
datMeta$Term_Accession_Number.1 = as.factor(datMeta$Term_Accession_Number.1)
datMeta$Study="PD.Riley"

#Load quantile normaliezed micorarray data
#Downlaod from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28475&format=file&file=GSE28475%5Fquantile%5Fnormalized%2Etxt%2Egz
datExpr = read.delim(file="./raw_data/Microarray/PD_Riley/GSE54282_series_matrix.txt",comment="!")
rownames(datExpr) = datExpr[,1]
datExpr=datExpr[,-1]
datExpr <- log(datExpr,base=2)
datExpr <- normalizeQuantiles(datExpr)


idx = match(colnames(datExpr), rownames(datMeta))
datMeta = datMeta[idx,]

#Initial QC Plots
boxplot(datExpr, range = 0, col= as.numeric(datMeta$disease_state), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");  legend("top", legend=levels(datMeta$disease_state),pch=19,col=c(1:2))
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$disease_state)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$disease_state), col=c(1:length(levels(datMeta$disease_state))), pch=16, cex=0.8)
plot(hclust(dist(t(datExpr)),method="average"))


##Annotate Probes
#ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
#identifier <- "entrezgene"
#getinfo <- c("affy_hugene_2_0_st_v1", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
#geneDat<- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
#geneDat<-getBM(attributes=c("affy_hugene_2_0_st_v1", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position"), filters="entrezgene", values=rownames(datExpr),
#               mart = ensembl)
geneDat <- read.delim("./raw_data/Microarray/PD_Riley/GPL17047_HuGene10stv1_Hs_ENTREZG_probe_tab.txt",header=T)
geneDat <- rownames(datExpr)
geneDat <- data.frame(geneDat)
names(geneDat) <- gsub("geneDat","Probe",names(geneDat))
geneDat$entrenzgene <- gsub("_at","",rownames(datExpr))

idx = match(rownames(datExpr), geneDat$Probe)
datProbes = cbind(rownames(datExpr), geneDat[idx,])

colnames(datExpr) = rownames(datMeta)
save(file="./working_data/Microarray/01_Normalized/Microarray_PD_Riley_normalized.RData", datMeta, datExpr, datProbes)
pdf("./results/figures/MicroarrayQC/PD_Riley_QC.pdf", width=11,height=8.5)
par(mfrow=c(3,4))
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$disease_state), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
legend("bottomleft",pch=16, legend = levels(datMeta$disease_state), col = 1:2)
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

model = model.matrix(~disease_state+age+Sex, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_PD_Riley_normalized_balanced.RData", datMeta, datExpr, datProbes,model)

#Post-Normalization QC
boxplot(datExpr, range = 0, col= as.numeric(datMeta$disease_state), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");   
legend("topright", legend=levels(datMeta$disease_state), col = 1:2, pch=19)
i = 1; plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease_state[i]), main="Hist of Log2 Exp", xlab = "log2 exp");   for(i in 2:dim(datExpr)[2]) {     lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$disease_state[i]),) } ;   legend("topright", levels(datMeta$disease_state), cex=0.7, text.col = 1:2)

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$disease_state)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$disease_state), col=as.numeric(datMeta$disease_state), pch=16, cex=0.8)

# Plot potential Colinearities
plot(datMeta$disease_state, ylim=c(0,20), ylab="Number", main="Subjects")
#A = anova(lm(as.numeric(datMeta$Term_Accession_Number.1) ~ datMeta$disease_state)); p = A$"Pr(>F)"[1];   plot(datMeta$disease_state ~ datMeta$datMeta$Term_Accession_Number.1, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$disease_state) ~ as.factor(datMeta$Term_Accession_Number.1))); p = A$"Pr(>F)"[1];   plot(datMeta$disease_state ~ as.factor(datMeta$Term_Accession_Number.1), main=paste("Chip, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
#A = anova(lm(as.numeric(datMeta$disease_state) ~ as.factor(datMeta$ChipPosition))); p = A$"Pr(>F)"[1];   plot(datMeta$disease_state ~ as.factor(datMeta$ChipPosition), main=paste("Chip Position, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$disease_state)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$disease_state, main=paste("Sex, p=1"), ylab="", xlab="")
A = anova(lm(as.numeric(datMeta$age) ~ datMeta$disease_state)); p = A$"Pr(>F)"[1];   plot(datMeta$age ~ datMeta$disease_state, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$disease_state)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$disease_state, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$disease_state)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$disease_state, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");

# Cluster Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$age),max(datMeta$age)))
plotDendroAndColors(tree, colors = cbind(as.numeric(datMeta$disease_state), as.numeric(datMeta$Term_Accession_Number.1), sex_col, age_col), disease_stateLabels = c("disease_state", "Batch","Sex", "Age"), cex.colorLabels=0.6, cex.dendroLabels=0.5) 

#PCA for sample
pca <-prcomp(t(datExpr),scale=F)
group <-factor(datMeta$disease_state)
colour_group <-cm.colors(length(unique(group)))
colour <-colour_group[as.numeric(factor(group))]
colour
group2 <-data.frame(group)
data2 <-cbind(t(datExpr),group2)
autoplot(pca,data=data2,colour='group')
dev.off()

# CollapseRows Probes --> Genes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
geneAnno<-getBM(attributes=c("ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position"), filters="entrezgene", values=datProbes$entrenzgene,
                mart = ensembl)
names(datProbes) <- gsub("entrenzgene","entrezgene",names(datProbes))
datProbes <- merge(geneAnno,datProbes,by="entrezgene")
datProbes <- datProbes[!duplicated(datProbes$Probe),]
rownames(datProbes) <- datProbes$`rownames(datExpr)`
realGenes = !is.na(datProbes$ensembl_gene_id)  
table(realGenes)
datExpr1 = datExpr[realGenes,]; datProbes1 = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$`rownames(datExpr)`) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$`rownames(datExpr)`)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Regress all Technical and Non-DX biological covariates
#datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
X = model.matrix(~disease_state+age+Sex, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,3:4]) %*% (as.matrix(beta[3:4,]))) 
datExpr = datExpr - t(to_regress)

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_PD_Riley_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)
