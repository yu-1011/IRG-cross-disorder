#######Calculation of module preservation
load("./working_data//NetworkAnalysis/eachdisordernet0802.rdata")
setLabels = c("disorders","background");
multiExpr = list(disorders = list(data =ad_datExpr),background = list(scz_datExpr,asd_datExpr,bd_datExpr,mdd_datExpr,pd_datExpr));
multiColor = list(disorders=ad_net$colors,background = ad_bg_net$colors)

#ref=AD
for (a in c("ad","pd","asd")){
for (i in c("scz","asd","bd","mdd","pd","ad")){
multiExpr = list(disorders = list(data=get(paste(a,"_datExpr",sep=""))), background = list(data=get(paste(i,"_datExpr",sep = ""))));
multiColor = list(disorders=get(paste(a,"_net",sep=""))$colors,background = get(paste(i,"_net",sep = ""))$colors)

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
#save(mp, file = "./processed_data/WGCNA/SCZ_modulePreservationin.RData");
#load("./processed_data/WGCNA/SCZ_modulePreservationin.RData")
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
z <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
z$disorder <- "SCZ"
write.table(z,paste("./results/tables/Manuscript/",i,"vs_",a,"_module_z-summary.txt",sep=""))
}}

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
pdf(file="./results/figures/WGCNA/mdd_modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();

# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
# Start the plot: open a suitably sized graphical window and set sectioning and margins.
pdf(file="./results/figures/WGCNA/mdd_modulePreservation.pdf", wi=10, h=5)
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
dev.off()

moduleColors <- scz_net$colors
immunegene <- read.table("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList",header=T,sep="\t")
table.p = matrix(NA, nrow=length(unique(moduleColors)), ncol=1)
rownames(table.p) = unique(moduleColors); 
colnames(table.p) = c("p.value")
table.or = table.p
hgnc = colnames(datExpr)
for (m in unique(moduleColors)) {
  for(e in colnames(table.p)) {
    #f = ORA(hgnc[CNS_moduleColors==m],cns.smr$Chinese.specific.eQTL.genes, hgnc, cns.smr$Chinese.specific.eQTL.genes)
    table.or[m,e] = length(intersect(hgnc[moduleColors==m],immunegene$ensembl))
    table.p[m,e] = dhyper(length(intersect(hgnc[moduleColors==m],immunegene$ensembl)),length(intersect(hgnc,immunegene$ensembl)),length(hgnc)-length(intersect(hgnc,immunegene$ensembl)),length(hgnc[moduleColors==m]))
    #table.or[m,e] = as.numeric(f[[1]])
    #table.p[m,e] = as.numeric(f[[2]])
  }}
table.p <- as.data.frame(table.p)
table.p$fdr = p.adjust(table.p$p.value,"bonferroni")
table.p$disorder = "scz"
write.table(table.p,"./results/tables/Manuscript/MDD_module_immuneEnrich.txt")
