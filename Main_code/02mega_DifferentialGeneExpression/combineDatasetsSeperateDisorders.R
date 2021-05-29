##3a_combineDatasets
## This script combines microarray studies together and performs ComBat to normalize them
## QC plots are made


rm(list=ls()); options(stringsAsFactors=F)

#source("http://bioconductor.org/biocLite.R")
library(ggplot2);library(sva); library(WGCNA)
rootdir = "/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/"
setwd(rootdir)

#par(mfrow=c(2,2))

files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_")
#files = files[grep("SCZ",files)]
files = files[!grepl("IBD",files)]
multiExpr = vector(mode="list",length = length(files))
for( i in 1:length(files)) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/",files[[i]],sep=""))
  print(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/",files[[i]],sep=""))
  multiExpr[[i]]$datExpr= datExpr
  multiExpr[[i]]$datMeta= datMeta
  multiExpr[[i]]$datProbes= datProbes
  rm(datExpr); rm(datMeta); rm(datProbes)
}

genes = rownames(multiExpr[[1]]$datExpr)
for(i in 2:length(multiExpr)) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))
idx = match(genes, rownames(multiExpr[[2]]$datProbes))
datProbes = multiExpr[[2]]$datProbes[idx, c("ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")]

all_datExpr = data.frame(row.names = genes)
all_datMeta = data.frame(matrix(NA, nrow=0, ncol=1));


for(i in 1:length(multiExpr)) {
  all_datExpr = cbind(all_datExpr, multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),])
  datMetaNA = as.data.frame(matrix(NA, ncol=10, nrow=nrow(multiExpr[[i]]$datMeta)))
  idx = pmatch(c("Study", "Subject", "Group", "Region", "Age", "Sex", "PMI", "pH", "RIN", "RNAdeg"),colnames(multiExpr[[i]]$datMeta))
  datMetaNA[,which(!is.na(idx))] = multiExpr[[i]]$datMeta[,na.omit(idx)]
  all_datMeta=rbind(all_datMeta, datMetaNA)
}

colnames(all_datMeta) = c("Study", "Subject", "Group", "Region", "Age", "Sex", "PMI", "pH", "RIN", "RNAdeg")
all_datMeta$RNAdeg[all_datMeta$Study=="SCZ.BD.Chen"]= NA #3' bias is not compatible from this study to the others (array platforM)
all_datMeta$Name <- colnames(all_datExpr)
all_datExpr <- all_datExpr[,!duplicated(colnames(all_datExpr))]
all_datMeta <- all_datMeta[!duplicated(all_datMeta$Name),]
rownames(all_datMeta) <- all_datMeta$Name
all_datMeta$Sex <- gsub("Female","F",all_datMeta$Sex)
all_datMeta$Sex <- gsub("female","F",all_datMeta$Sex)
all_datMeta$Sex <- gsub("Male","M",all_datMeta$Sex)
all_datMeta$Sex <- gsub("male","M",all_datMeta$Sex)
all_datMeta$Sex <- as.factor(all_datMeta$Sex)
all_datMeta$Group <- gsub("control","CTL",all_datMeta$Group)
all_datMeta$Group <- gsub("Control","CTL",all_datMeta$Group)
all_datMeta$Group <- gsub("normal","CTL",all_datMeta$Group)
all_datMeta$Group <- gsub("Parkinson's_disease","PD",all_datMeta$Group)
all_datMeta$Group <- gsub("Affected","AD",all_datMeta$Group)
all_datMeta$Group <- gsub("Alzheimer's_Disease","AD",all_datMeta$Group)
all_datMeta$Group <- as.factor(all_datMeta$Group)

#datMeta = subset(all_datMeta,Group=="SCZ"|Group=="BD"|Group=="CTL")
#datMeta$Group <- factor(datMeta$Group,levels =c("SCZ","BD","CTL"))
#datExpr = all_datExpr[,rownames(datMeta)]

datMeta = all_datMeta
datExpr = all_datExpr


##QC Pre-Combat
pdf("./results/figures/MicroarrayQC/CombinedDatasetPre-Combat.pdf")
datMeta$Study = as.factor(datMeta$Study)
sex_col = rep("blue", times = nrow(datMeta))
sex_col[datMeta$Sex=="F"] = "pink"
age_col = numbers2colors(datMeta$Age, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Age, na.rm=T),max(datMeta$Age, na.rm=T)))
#ph_col = numbers2colors(datMeta$pH, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$pH, na.rm=T),max(datMeta$pH, na.rm=T)))
#pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI, na.rm=T),max(datMeta$PMI, na.rm=T)))
#rin_col = numbers2colors(datMeta$RIN, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RIN, na.rm=T),max(datMeta$RIN, na.rm=T)))
#rna_col = numbers2colors(datMeta$RNAdeg, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$RNAdeg, na.rm=T),max(datMeta$RNAdeg, na.rm=T)))

plot(density(datExpr[,1]), xlim=c(-5,20), ylim=c(0, 0.5), col = as.numeric(datMeta$Study[1]), xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Pre-Combat")
for(i in 2:dim(datExpr)[[2]])
  lines(density(datExpr[,i]), xlim=c(0,20), col = as.numeric(datMeta$Study[i]))  
legend("topright", (levels(datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Study)), pch=20, main="Multidimensional Scaling Plot\nPre-ComBat", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topleft", (levels(datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

tree = hclust(dist(t(datExpr)), method = "average")
par(mfrow=c(1,1))
plotDendroAndColors(tree, cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Study), sex_col, age_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age"), cex.colorLabels=0.8, cex.dendroLabels=0.15,
                    main="Dendrogram\nPre-Combat")
dev.off()


#Normalize by Study
mod = model.matrix(~Sex+Age+Group, data=datMeta)
batch = as.factor(datMeta$Study)
datExpr = ComBat(datExpr, batch=batch, mod=mod)

##QC - PostCombat
pdf("./results/figures/MicroarrayQC/CombinedDatasetPost-Combat.pdf")
par(mfrow=c(2,2))
plot(density(datExpr[,1]), col = as.numeric(datMeta$Study[1]), , xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Post-Combat")
for(i in 2:dim(datExpr)[[2]])
  lines(density(datExpr[,i]), xlim=c(0,16), ylim=c(0,0.3), col = as.numeric(datMeta$Study[i]))  
legend("topright", levels(datMeta$Study), col=c(1:8), pch=16,cex=0.7)

#MDS Plot
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Study)), pch=20, main="MDS: Study", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=as.numeric(as.factor(datMeta$Group)), pch=16, main="MDS: Group",  xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
legend("bottomleft", levels(datMeta$Group), col=c(1:length(levels(datMeta$Group))), pch=16, cex=0.8)
plot(mds$points, col=sex_col, pch=16, main="MDS - Sex", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=age_col, pch=16, main="MDS - Age", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=pmi_col, pch=16, main="MDS - PMI", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=ph_col, pch=16, main="MDS - pH", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=rin_col, pch=16, main="MDS - RIN", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
#plot(mds$points, col=rna_col, pch=16, main="MDS - RNA", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 

#Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
plotDendroAndColors(tree, cbind(as.numeric(datMeta$Group), as.numeric(datMeta$Study), sex_col, age_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age", "pH", "PMI", "RIN", "RNA"), cex.colorLabels=0.8, cex.dendroLabels=0.15,
                    main="Dendrogram\nPost-Combat")

dev.off()

save(file = "./working_data//NetworkAnalysis/All_datasets_combined_190410.RData", datExpr,datMeta,multiExpr)
write.table(datExpr,file="./working_data/SeparateCombinedData/ASD_Expr.txt",quote = F)
write.table(datMeta,file="./working_data/SeparateCombinedData/ASD_Meta.txt",quote = F)
