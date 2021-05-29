##1b_Microarray_AD_Stephan.R

rm(list=ls()); options(stringsAsFactors = F)

library(limma); library(WGCNA); library(biomaRt); library(sva)
setwd("/Users/normacy/Desktop/2018immuneGene/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master")

#Load & clean metaData
datMeta = read.delim("./raw_data/Microarray/AD_GSE63063_Jamie_blood//E-GEOD-63063.sdrf.txt")
rownames(datMeta) = datMeta$Comment..Sample_title.
datMeta$Subjects = datMeta$Comment..Sample_title.

datMeta$Disease_State = gsub("MCI","AD",datMeta$Characteristics..status.)
datMeta$Disease_State = factor(datMeta$Disease_State,levels =  c("CTL", "AD"))

datMeta$Sex = as.factor(datMeta$Characteristics..sex.)
datMeta$Age = as.numeric(datMeta$Characteristics..age.)
datMeta$Ethnicity = as.factor(datMeta$Characteristics..ethnicity.)
#datMeta$Cell_Type = as.factor(datMeta$Cell_Type)
datMeta$Study="AD.Jamie"
datMeta_all = datMeta
  
#Load quantile normaliezed micorarray data
#Downlaod from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28475&format=file&file=GSE28475%5Fquantile%5Fnormalized%2Etxt%2Egz
datExpr = read.delim(file="./raw_data/Microarray/AD_GSE63061_Jamie_blood/GSE63061_normalized (1).txt",comment="!",header=F)
datExpr<-datExpr[which(duplicated(datExpr$V1)!=T),]
rownames(datExpr) = datExpr[,1]
datExpr=datExpr[,-1]
colnames(datExpr) = datExpr[1,]
datExpr=datExpr[-1,]
datExpr <- log(datExpr,base=2)
datExpr <- normalizeQuantiles(datExpr)
X <- apply(datExpr, 2, as.numeric)
rownames(X) <- rownames(datExpr) 
datExpr <- X

idx = match(colnames(datExpr), datMeta_all$Subjects)
datMeta = datMeta_all[idx,]

#Initial QC Plots
boxplot(datExpr,col= as.numeric(datMeta$Disease_State), xaxt='n', xlab = "Array", main ="Boxplot", ylab = "Intensity");legend("top",legend=levels(datMeta$Disease_State),pch=19,col=1:2);legend("top",legend=levels(datMeta$Disease_State),col=1:2)
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Disease_State)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Disease_State), col=c(1:length(levels(datMeta$Disease_State))), pch=16, cex=0.8)
plot(hclust(dist(t(datExpr)),method="average"))


##Annotate Probes
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
identifier <- "illumina_humanht_12_v4"
getinfo <- c("illumina_humanht_12_v4", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
#geneDat<- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
geneDat<-getBM(attributes=c("illumina_humanht_12_v4", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position"), filters="illumina_humanht_12_v4", values=rownames(datExpr),
                          mart = ensembl)

idx = match(rownames(datExpr), geneDat$illumina_humanht_12_v4)
datProbes = cbind(rownames(datExpr), geneDat[idx,])

colnames(datExpr) = rownames(datMeta)
save(file="./working_data/Microarray/01_Normalized/Microarray_AD_Jamie_2_blood_normalized.RData", datMeta, datExpr, datProbes)
pdf("./results/figures/MicroarrayQC/AD_Jamie_2_QC.pdf", width=11,height=8.5)
par(mfrow=c(3,4))
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Disease_State), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
legend("bottomleft",pch=16, legend = levels(datMeta$Disease_State), col = 1:2)
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

model = model.matrix(~Disease_State+Age+Sex+Ethnicity, data=datMeta)
save(file="./working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_AD_Jamie_2_blood_normalized_balanced.RData", datMeta, datExpr, datProbes,model)

#Post-Normalization QC
boxplot(datExpr, range = 0, col= as.numeric(datMeta$Disease_State), xaxt='n', xlab = "Array", main ="Boxplot Post-Normalization", ylab = "Intensity");legend("top",legend=levels(datMeta$Disease_State),pch=19,col=1:2)   
#legend("topright", legend=levels(datMeta$Disease_State), col = 1:2, pch=19)
i = 1
plot(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Disease_State[i]), main="Hist of Log2 Exp", xlab = "log2 exp"); 
for(i in 2:dim(datExpr)[2]) {lines(density((datExpr[,i]), na.rm=T), col = as.numeric(datMeta$Disease_State[i]),)} 
legend("topright", levels(datMeta$Disease_State), cex=0.7,text.col = 1:2)

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Disease_State)), pch=16, main="MDS Plot", asp=1, xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep=""));   legend("bottomright", levels(datMeta$Disease_State), col=as.numeric(datMeta$Disease_State), pch=16, cex=0.8)

# Plot potential Colinearities
plot(datMeta$Disease_State, ylim=c(0,100), ylab="Number", main="Subjects")
#A = anova(lm(as.numeric(datMeta$Age) ~ datMeta$Disease_State)); p = A$"Pr(>F)"[1];   plot(datMeta$Disease_State ~ datMeta$Batch, main=paste("Batch, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$Disease_State) ~ as.factor(datMeta$Organ_Region))); p = A$"Pr(>F)"[1];   plot(datMeta$Disease_State ~ as.factor(datMeta$Organ_Region), main=paste("Organ_Region, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$Disease_State) ~ as.factor(datMeta$Ethnicity))); p = A$"Pr(>F)"[1];   plot(datMeta$Disease_State ~ as.factor(datMeta$Ethnicity), main=paste("Ethnicity, p=", signif(p,2)), ylab="", xlab="", xlim=c(0,1))
A = anova(lm(as.numeric(datMeta$Sex) ~ datMeta$Disease_State)); p = A$"Pr(>F)"[1];   plot(datMeta$Sex ~ datMeta$Disease_State, main=paste("Sex, p=1"), ylab="", xlab="")
A = anova(lm((datMeta$Age) ~ datMeta$Disease_State)); p = A$"Pr(>F)"[1];   plot(datMeta$Age ~ datMeta$Disease_State, main=paste("Age, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$PMI) ~ datMeta$Disease_State)); p = A$"Pr(>F)"[1];   plot(datMeta$PMI ~ datMeta$Disease_State, main=paste("PMI, p=", signif(p,2)), ylab="", xlab="")
#A = anova(lm(as.numeric(datMeta$RIN) ~ datMeta$Disease_State)); p = A$"Pr(>F)"[1];   plot(datMeta$RIN ~ datMeta$Disease_State, main=paste("RIN, p=", signif(p,2)), ylab="", xlab="");

# Cluster Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
sex_col = as.numeric((datMeta$Sex)); sex_col[sex_col==2] = "blue"; sex_col[sex_col==1]="pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age),max(datMeta$Age)))
plotDendroAndColors(tree, colors = cbind(as.numeric(datMeta$Disease_State), sex_col, age_col), Disease_StateLabels = c("Disease_State","Sex","Age"), cex.colorLabels=0.6, cex.dendroLabels=0.5)
dev.off()



# CollapseRows Probes --> Genes
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$illumina_humanht_12_v4) 
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$illumina_humanht_12_v4)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Regress all Technical and Non-DX biological covariates
#datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
X = model.matrix(~Disease_State, data=datMeta)
datExpr <- Y[,which(rownames(X) %in% colnames(Y))]
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,1:2]) %*% (as.matrix(beta[1:2,]))) 
datExpr = Y - t(to_regress)

#arrangedatMeta
datMeta$Group <- datMeta$Disease_State

save(file = "./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AD_Jamie_2_blood_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)
