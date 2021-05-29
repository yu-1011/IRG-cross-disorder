#3c_networkDendroAndParams.R

rm(list=ls()); options(stringsAsFactors = F)

library(WGCNA)
#Load combined datasets
load("./working_data/NetworkAnalysis/All_datasets_combined_081518_11245x625.RData")

#load TOM comptued by rWGCNA
load("./working_data/NetworkAnalysis/Brain-combined-consensusTOM_final0815.RData")
geneTree = hclust(1-as.dist(tom), method="average")


# Iterate WGCNA parameters for robustness -- this takes a while
if(TRUE) {
   colors = vector(mode="list")
  labels = vector(mode="list")
  for (pam in c(FALSE,TRUE)) {
    for (minModSize in c(50,100, 200)) {
      for (dthresh in c(0.1, 0.2)) {
        for(ds in c(0:4)) { 
            print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
            
            tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-as.dist(tom)))
            merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=dthresh)
            colors = cbind(colors, labels2colors(merged$colors))
            
            labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        }
      }
    }
  }
  
  plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE, dendroLabels=FALSE, main="Dendrogram", cex.colorLabels=0.5)
  save(file="./working_data/NetworkAnalysis//WGCNA_diffParams0815.rda", geneTree, colors, labels)
  
  colors2 = colors
  colors2[,seq(1:30)] = colors[,seq(1,60,by=2)]
  colors2[,seq(31:60)] = colors[,seq(2,60,by=2)]
  
  plotDendroAndColors(geneTree,colors,addGuide=T,dendroLabels=F)
  
  for (i in 1:60){
    ci = as.character(colors[,i])
    c_new = matchLabels(ci, c_ref)
    colors[,i] = c_new
  }
  
  colors = cbind(colors[,15], colors)
  labels = c("Final Modules", labels)
  
  pdf("./results/figures/WGCNA/WGCNA_diffParams0815.pdf",width=6,height=8)
  plotDendroAndColors(geneTree,colors,groupLabels = labels,addGuide=T,dendroLabels=F,cex.colorLabels=0.3)
  dev.off()
  
}



# Finalized Parameters
# --------------------
# Parameters to Use: "DS=4,MMS=50,DCOR=0.2,PAM=FALSE"
wgcna_parameters = list(powers =  8)
wgcna_parameters$minModSize = 50
wgcna_parameters$minHeight = 0.2
wgcna_parameters$bsize = 26000  ##block size needs to be larger than dim(datExpr)[1]
wgcna_parameters$ds = 4  ##deep split parameter contorls number of modules
wgcna_parameters$networkType = "signed"    ## using signed networks
wgcna_parameters$corFnc = "bicor"
wgcna_parameters$pamStage = F

tree = cutreeHybrid(dendro = geneTree, minClusterSize= wgcna_parameters$minModSize, pamStage=wgcna_parameters$pamStage, cutHeight = 0.999, 
                    deepSplit=wgcna_parameters$ds, distM=as.matrix(1-as.dist(tom)))

merged = mergeCloseModules(exprData= t(datExpr), colors = tree$labels, cutHeight=wgcna_parameters$minHeight)

#merged$colors <- gsub("13","12",merged$colors)
#merged$colors <- gsub("14","13",merged$colors)
merged$colors <- as.numeric(merged$colors)

moduleColors = labels2colors(merged$colors)
colors = moduleColors
table(moduleColors)
length(table(moduleColors))

plotDendroAndColors(geneTree,moduleColors,groupLabels = "mod",cex.colorLabels = 0.5,addGuide=T,dendroLabels=F)


MEs = moduleEigengenes(expr = t(datExpr), moduleColors, softPower = wgcna_parameters$powers)

colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
save(MEs, moduleLabels, moduleColors, geneTree, file = "./working_data/NetworkAnalysis/Brain-Combined-signed-afterDP-networkConstruction-stepByStep0903.RData")

kMEtable = signedKME(t(datExpr),MEs$eigengenes)
tableS1 = data.frame(kMEtable[,paste0("kME", unique(moduleColors))])
colnames(tableS1) = paste0("kME.CD", 1:11, ".", unique(moduleColors))
tableS1 = cbind( data.frame(Module.Color=colors, Module.name = paste0("CD",merged$colors)), tableS1)
igenes = intersect(rownames(tableS1),rownames(immunegene))
tableS2 = tableS1[igenes,]
write.csv(file="./results/tables/Manuscript/TableS1 - brainCombined - kME table0903.csv", tableS1)
write.csv(file="./results/tables/Manuscript/TableS2 - brainCombinedimmuneGene - kME table0903.csv", tableS1)

save(file="./working_data/NetworkAnalysis/brainCombinedfinalizedNetwork_0903.RData", datExpr, datMeta, multiExpr, geneTree, colors, wgcna_parameters, colors,MEs, kMEtable)
