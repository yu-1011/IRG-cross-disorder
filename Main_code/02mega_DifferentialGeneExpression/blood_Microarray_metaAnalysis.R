##4a_BloodMicroarray_meta-Analysis
rm(list=ls()); options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R")

library(nlme);library(limma)
setwd("/Users/normacy/Desktop/immuneGene/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap-master")
immunegene <- read.table("/Users/normacy/Desktop/immuneGene/01data/ImmuneGeneList",header=T,row.names = 1)


##AD
load("./working_data/Microarray/blood/blood/Microarray_AD_Giovana_blood_normalized_CR_cleaned.RData")
AD_datExpr = datExpr; rm(datExpr)
AD_datMeta = datMeta; rm(datMeta)
AD_datProbes = datProbes; rm(datProbes)

genes= rownames(AD_datExpr)
AD_datExpr = AD_datExpr[match(genes,rownames(AD_datExpr)),]
AD_Group = factor(AD_datMeta$disease, levels=c("CTL", "AD"))
AD_Study = as.factor(AD_datMeta$Study)
AD_Subject = as.factor(AD_datMeta$Source.Name)
AD_datMeta <- data.frame(Study=AD_Study, Subject=AD_Subject, Group=AD_Group)
AD_meta = matrix(NA, nrow=nrow(AD_datExpr), ncol=3)
for(i in 1:nrow(AD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AD_datExpr[i,])
  tryCatch({
    AD_meta[i,] = summary(lme(expr~ Group,data = AD_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
AD_meta=as.data.frame(AD_meta)
colnames(AD_meta) = c("beta", "SE", "p")
rownames(AD_meta) = genes
AD_meta$fdr = p.adjust(AD_meta$p, "fdr")
AD_meta$symbol=AD_datProbes$external_gene_id[match(genes, AD_datProbes$ensembl_gene_id)]
#AD_meta= AD_meta[!apply(is.na(AD_meta),1,any),]
write.csv(file="./results/tables/Microarray_AD_metaanalysis_blood_082318.csv",AD_meta)

##SCZ
load("./working_data/Microarray/blood/blood/Microarray_SCZ_Horiuchin_blood_normalized_CR_cleaned.RData")
SCZ_datExpr = datExpr; rm(datExpr)
SCZ_datMeta = datMeta; rm(datMeta)
SCZ_datProbes = datProbes; rm(datProbes)

genes= rownames(SCZ_datExpr)
SCZ_datExpr = SCZ_datExpr[match(genes,rownames(SCZ_datExpr)),]
SCZ_Group = factor(SCZ_datMeta$Group, levels=c("CTL", "SCZ"))
SCZ_Study = as.factor(SCZ_datMeta$Study)
SCZ_Subject = as.factor(SCZ_datMeta$Sample_geo_accession)
SCZ_datMeta <- data.frame(Study=SCZ_Study, Subject=SCZ_Subject, Group=SCZ_Group)


SCZ_meta = matrix(NA, nrow=nrow(SCZ_datExpr), ncol=3)
for(i in 1:nrow(SCZ_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(SCZ_datExpr[i,])
  tryCatch({
    SCZ_meta[i,] = summary(lme(expr~ Group,data = SCZ_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
SCZ_meta=as.data.frame(SCZ_meta)
colnames(SCZ_meta) = c("beta", "SE", "p")
rownames(SCZ_meta) = genes
SCZ_meta$fdr = p.adjust(SCZ_meta$p, "fdr")
SCZ_meta$symbol=SCZ_datProbes$external_gene_id[match(genes, SCZ_datProbes$ensembl_gene_id)]
#SCZ_meta= SCZ_meta[!apply(is.na(SCZ_meta),1,any),]
write.csv(file="./results/tables/Microarray_SCZ_metaanalysis_blood_082318.csv",SCZ_meta)

##BD
load("./working_data/Microarray/blood/blood/Microarray_BD_Clelland_blood_normalized_CR_cleaned.RData")
BD_datExpr = datExpr; rm(datExpr)
BD_datMeta = datMeta; rm(datMeta)
BD_datProbes = datProbes; rm(datProbes)

genes= rownames(BD_datExpr)
BD_datExpr = BD_datExpr[match(genes,rownames(BD_datExpr)),]
BD_Group = factor(BD_datMeta$disease, levels=c("CTL", "BD"))
BD_Study = as.factor(BD_datMeta$Study)
BD_Subject = as.factor(BD_datMeta$Source.Name)
BD_datMeta <- data.frame(Study=BD_Study, Subject=BD_Subject, Group=BD_Group)


BD_meta = matrix(NA, nrow=nrow(BD_datExpr), ncol=3)
for(i in 1:nrow(BD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(BD_datExpr[i,])
  tryCatch({
    BD_meta[i,] = summary(lme(expr~ Group,data = BD_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
BD_meta=as.data.frame(BD_meta)
colnames(BD_meta) = c("beta", "SE", "p")
rownames(BD_meta) = genes
BD_meta$fdr = p.adjust(BD_meta$p, "fdr")
BD_meta$symbol=BD_datProbes$external_gene_id[match(genes, BD_datProbes$ensembl_gene_id)]
#BD_meta= BD_meta[!apply(is.na(BD_meta),1,any),]
write.csv(file="./results/tables/Microarray_BD_metaanalysis_blood_082318.csv",BD_meta)

#PD
load("./working_data/Microarray/blood/blood/Microarray_PD_Paolan_blood_normalized_CR_cleaned.RData")
PD_datExpr = datExpr; rm(datExpr)
PD_datMeta = datMeta; rm(datMeta)
PD_datProbes = datProbes; rm(datProbes)

genes= rownames(PD_datExpr)
PD_datExpr = PD_datExpr[match(genes,rownames(PD_datExpr)),]
PD_Group = factor(PD_datMeta$disease, levels=c("CTL", "PD"))
PD_Study = as.factor(PD_datMeta$Study)
PD_Subject = as.factor(PD_datMeta$Source.Name)
PD_datMeta <- data.frame(Study=PD_Study, Subject=PD_Subject, Group=PD_Group)


PD_meta = matrix(NA, nrow=nrow(PD_datExpr), ncol=3)
for(i in 1:nrow(PD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(PD_datExpr[i,])
  tryCatch({
    PD_meta[i,] = summary(lme(expr~ Group,data = PD_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
PD_meta=as.data.frame(PD_meta)
colnames(PD_meta) = c("beta", "SE", "p")
rownames(PD_meta) = genes
PD_meta$fdr = p.adjust(PD_meta$p, "fdr")
PD_meta$symbol=PD_datProbes$external_gene_id[match(genes, PD_datProbes$ensembl_gene_id)]
#PD_meta= PD_meta[!apply(is.na(PD_meta),1,any),]
write.csv(file="./results/tables/Microarray_PD_metaanalysis_blood_082318.csv",PD_meta)

load("./working_data/Microarray/blood/blood/Microarray_ASD_Kuwanon_blood_normalized_CR_cleaned.RData")
ASD_datExpr = datExpr; rm(datExpr)
ASD_datMeta = datMeta; rm(datMeta)
ASD_datProbes = datProbes; rm(datProbes)

genes= rownames(ASD_datExpr)
ASD_datExpr = ASD_datExpr[match(genes,rownames(ASD_datExpr)),]
ASD_Group = factor(ASD_datMeta$disease, levels=c("CTL", "ASD"))
ASD_Study = as.factor("ASD_blood")
ASD_Subject = as.factor(ASD_datMeta$Source.Name)
ASD_datMeta <- data.frame(Study=ASD_Study, Subject=ASD_Subject, Group=ASD_Group)


ASD_meta = matrix(NA, nrow=nrow(ASD_datExpr), ncol=3)
for(i in 1:nrow(ASD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(ASD_datExpr[i,])
  tryCatch({
    ASD_meta[i,] = summary(lme(expr~ Group,data = ASD_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
ASD_meta=as.data.frame(ASD_meta)
colnames(ASD_meta) = c("beta", "SE", "p")
rownames(ASD_meta) = genes
ASD_meta$fdr = p.adjust(ASD_meta$p, "fdr")
ASD_meta$symbol=ASD_datProbes$external_gene_id[match(genes, ASD_datProbes$ensembl_gene_id)]
#ASD_meta= ASD_meta[!apply(is.na(ASD_meta),1,any),]
write.csv(file="./results/tables/Microarray_ASD_metaanalysis_blood_082318.csv",ASD_meta)


