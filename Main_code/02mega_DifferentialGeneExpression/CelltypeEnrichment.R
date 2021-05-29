rm(list=ls())
rootdir = "/Users/normacy/Desktop/immuneGene/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master"
setwd(rootdir)
library(pSI); library(gProfileR); library(gplots); library(biomaRt); library(WGCNA)

load("./working_data/NetworkAnalysis/All_datasets_combined_081518_11245x625.RData")
load("./working_data/NetworkAnalysis/brainCombinedfinalizedNetwork_0815.RData")


eigmat = MEs$eigengenes; colnames(eigmat) = gsub("ME","",colnames(eigmat))
kME = signedKME(t(datExpr), MEs$eigengenes); colnames(kME) = gsub("kME", "", colnames(kME))


##Download human transcriptome data from:
##Zhang, Y. et al. Purification and Characterization of Progenitor and Mature Human Astrocytes Reveals Transcriptional and Functional Differences with Mouse. Neuron 89, 37–53 (2016).
##Supplemental Table 3, "Human data only" sheet
##---
if(FALSE) {
  zhang.datExpr = read.csv("./raw_data/Annotations/datExpr.zhangHuman.avgForPSI.csv",skip=3,nrow=23223,head=F) 
  zhang.datMeta = data.frame(row.names=1:41,t(read.csv("./raw_data/Annotations/nn.4063-S12.csv",nrow=3,head=F)[,-1]))
  zhang.datMeta$CellType = NA
  zhang.datProbes = data.frame(symbol=zhang.datExpr$V1)
  zhang.datExpr = zhang.datExpr[,-1]
  colnames(zhang.datMeta) = c("X1", "Age", "Gender","CellType")
  zhang.datMeta$CellType[15:26] = "Astrocyte"
  zhang.datMeta$CellType[27] = "Neuron"
  zhang.datMeta$CellType[28:32] = "Oligo"
  zhang.datMeta$CellType[33:35] = "Microglia"
  zhang.datMeta$CellType[36:37] = "Endothelial"
  zhang.datMeta$CellType[38:41] = "WholeCortex"
  
  zhang.datExpr2 = data.frame(matrix(NA, nrow=nrow(zhang.datExpr), ncol=5)); colnames(zhang.datExpr2)=  c("Neuron", "Astrocyte", "Oligo", "Microglia","Endothelial")
  zhang.datExpr2$Neuron = zhang.datExpr[,which(zhang.datMeta$CellType=="Neuron")]
  for(cell in colnames(zhang.datExpr2)[2:5]) {
    zhang.datExpr2[,cell] = apply(zhang.datExpr[,which(zhang.datMeta$CellType==cell)],1,mean)  
  }
  
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
  bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id"), filters = "hgnc_symbol", values = zhang.datProbes$symbol, mart=ensembl)
  
  zhang.datProbes = data.frame(symbol=zhang.datProbes$symbol, ensg=bm$ensembl_gene_id[match(zhang.datProbes$symbol, bm$external_gene_id)])
  cr =collapseRows(zhang.datExpr2, rowGroup = zhang.datProbes$ensg, rowID=1:nrow(zhang.datExpr2))
  zhang.datExpr = cr$datETcollapsed
  rm(zhang.datExpr2, zhang.datProbes)
  save(file=".")
}


# Calculate Cell-Type specificity of modules
#Zhang using pSI
zhang.datExpr = read.csv("./raw_data/Annotations//datExpr.zhangHuman.avgForPSI.csv",row.names=1)[,-1]
set.seed(100)
pSI.output = specificity.index(pSI.in=zhang.datExpr,bts=100,p_max=.1, e_min=0.3); 
pSI.count(pSI.output)

cell.p.zhang = matrix(NA, 11,5);  rownames(cell.p.zhang) = unique(colors)
colnames(cell.p.zhang) = colnames(pSI.output)


for(mod in rownames(cell.p.zhang)) {
  f = fisher.iteration(pSI.output, rownames(datExpr)[colors==mod],p.adjust = F)
  cell.p.zhang[mod,] = f$`0.05 - nominal`
}

cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
to_plot = cell.p.zhang.fdr[c("greenyellow", "green", "turquoise", "blue"),]
to_plot = cell.p.zhang.fdr

dendro.col = as.dendrogram(hclust(as.dist(1-bicor(zhang.datExpr)), method="average"))
#denro.row= as.dendrogram(hclust(as.dist(1-bicor(eigmat[,c("greenyellow", "green", "turquoise", "blue")])),method="average"))
denro.row= as.dendrogram(hclust(as.dist(1-bicor(eigmat)),method="average"))


pdf("./results/figures/Manuscript/Fig3F-CellType.pdf",width=6,height=5)
heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000],
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none",
          colsep=0:7,rowsep=0:8,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5,
          Rowv=denro.row, Colv=dendro.col,
          key=T,key.xlab="-log10(P)", cellnote=signif(to_plot,1), notecex=.8, notecol="black",main="Enrichment")
dev.off()


out = cell.p.zhang.fdr
colnames(out) = paste("Enrichment.", colnames(out), ".FDR",sep="")
write.csv(out,file="./results/tables/Manuscript/TableS2-CellType.csv")

##
library(gProfileR); library(gplots); library(biomaRt); library(WGCNA);library(gplots)

cell.rujia = matrix(NA, length(unique(colors)),4);  rownames(cell.rujia) = unique(colors)
colnames(cell.rujia) = c("astrocyte","microglia","neuron","oligodendrocyte")
moduleInfo = data.frame(rownames(datExpr),colors)
colnames(moduleInfo) = c("ensembl_gene_id","colors")
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org") 
bm = getBM(attributes = c("ensembl_gene_id", "external_gene_id"), filters = "ensembl_gene_id", values = moduleInfo$rownames.datExpr., mart=ensembl)
a <- merge(moduleInfo,bm)

celltype <- read.csv("/Users/normacy/Downloads/collectedMarker.csv",header=T)
astrocyte <- subset(celltype,cell_type1=="astrocyte")
microglia <- subset(celltype,cell_type1=="microglia")
neuron <- subset(celltype,cell_type1=="neuron")
oligodendrocyte <- subset(celltype,cell_type1=="oligodendrocyte")
immunegeneFunction <- read.table("/Users/normacy/Desktop/immuneGene/01data/ImmuneGeneFunction",header=T,sep="\t")
source("./code/00_scripts/fisher_overlap.R")



for(mod in rownames(cell.rujia)){
  for (i in colnames(cell.rujia)){
  cell.rujia[mod,i] = ORA(subset(celltype,cell_type1==i)$gene.1, a$external_gene_id[colors==mod],celltype$gene.1,a$external_gene_id)[[2]]
  }
}

for(mod in rownames(cell.rujia)){
  for (i in colnames(cell.rujia)){
    cell.rujia[mod,i] = ORA(subset(celltype,cell_type1==i)$gene.1,a$external_gene_id[colors==mod],celltype$gene.1,multiExpr[[1]]$datProbes$`Associated Gene Name`)[[2]]
  }
}

for(mod in rownames(cell.rujia)) {
  for (i in colnames(cell.rujia)){
  o = intersect(subset(celltype,cell_type1==i)$gene.1,a$external_gene_id[colors==mod])
  m = subset(celltype,cell_type1==i)$gene.1
  n = length(a$external_gene_id)-length(subset(celltype,cell_type1==i)$gene.1)
  k = a$external_gene_id[colors==mod]
  cell.rujia[mod,i] = phyper(length(o),length(m),n,length(k))
  }
}

to_plot = as.numeric(cell.rujia)
dim(to_plot) =  dim(cell.rujia)
rownames(to_plot) = rownames(cell.rujia)
colnames(to_plot) = colnames(cell.rujia)

heatmap.2(-log10(to_plot),col=blueWhiteRed(1000,1)[500:1000],
          scale="none",trace="none",cexRow = 0.8,cexCol = .8, density.info = "none",
          colsep=0:7,rowsep=0:8,sepcolor="grey",sepwidth=c(0.02,0.02),
          srtCol=45,offsetRow=0,offsetCol=-0.5,
          key=T,key.xlab="-log10(P)", cellnote=signif(to_plot,1), notecex=.8, notecol="black",main="Enrichment")
dev.off()

