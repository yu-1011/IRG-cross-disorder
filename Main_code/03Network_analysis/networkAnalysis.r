#3b_networkAnalysis.R

rm(list=ls())
library(WGCNA)
rootdir = "/Users/normacy/Desktop/immuneGene/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master"
setwd(rootdir)
load("./working_data/NetworkAnalysis/All_datasets_combined_081518_11245x625.RData")

#######construct WGCNA unsigned network stepByStep###########

#power selection
## ----------------
multiExpr = vector(mode="list", length=1)
datExpr = as.data.frame(t(datExpr))
bsize = 1000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

#unsigned
sft<-pickSoftThreshold(datExpr,powerVector = powers, corFnc = bicor, networkType = "signed", verbose = 5)
pdf(file = "./results/figures/WGCNA/signed_threshold.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Co-expression similarity and adjacency
adj<-adjacency(datExpr=datExpr,type='unsigned',power=3,corFnc='bicor')

# Turn adjacency into topological overlap
tom<-TOMsimilarity(adj, TOMType = "signed", verbose = 1)
save(adj,tom,file='./working_data/NetworkAnalysis/Brain-combined-consensusTOM_final190410.RData')

# Call the hierarchical clustering function
dissTOM = 1-tom
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
pdf(file="/zs32_2/home/ychen/2017GWAS/CNSBraingenetree0607.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2,pamStage = F,pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file="./results/figures/WGCNA/UnsignedNetwork0410.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file="/zs32_2/home/ychen/2017GWAS/CNSBrainUnsignedMetreeimmune.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr,dynamicColors, relabel = T,cutHeight = MEDissThres, verbose = 4)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
pdf(file="/zs32_2/home/ychen/2017GWAS/CNSBrainUnsignedgeneDendro0607.pdf",wi=9,he=6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "./working_data/NetworkAnalysis/Brain-Combined-signed-networkConstruction-immunegene0411.RData")

#Save module info
#输出module信息
table<-moduleLabels
informatio<-data.frame(gene=rownames(t(datExpr)),module=table,color=labels2colors(moduleColors))
dim(informatio)
write.csv(informatio,file="./results/tables/Brain-Combined-signed-networkConstruction-geneinfo0815.csv")

#首先针对所有基因画热图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 8); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")



