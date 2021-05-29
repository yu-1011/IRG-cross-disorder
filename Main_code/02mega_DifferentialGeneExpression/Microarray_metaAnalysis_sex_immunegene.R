##2a_Microarray_meta-Analysis
rm(list=ls()); options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R")

library(nlme)
setwd("/Users/normacy/Desktop/2018immuneGene/ImmuneGeneAnalysis/")
immunegene <- read.table("/Users/normacy/Desktop/2018immuneGene/01data/ImmuneGeneList",header=T,row.names = 1)

#AD
n=2
AD_multiExpr = vector(mode="list",length = 2)
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned//Microarray_AD_Blalock_normalized_CR_cleaned.Rdata")
i=1; AD_multiExpr[[i]]$datExpr = datExpr; AD_multiExpr[[i]]$datMeta = datMeta ;AD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AD_Stephan_normalized_CR_cleaned.RData")
i=2; AD_multiExpr[[i]]$datExpr = datExpr; AD_multiExpr[[i]]$datMeta = datMeta ;AD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)

genes= intersect(intersect(rownames(AD_multiExpr[[1]]$datExpr), rownames(AD_multiExpr[[2]]$datExpr)),rownames(immunegene))

for(i in 1:n) AD_multiExpr[[i]]$datExpr = AD_multiExpr[[i]]$datExpr[match(genes,rownames(AD_multiExpr[[i]]$datExpr)),]

AD_datExpr = data.frame(row.names = genes)
AD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)

i=1;
AD_datExpr = cbind(AD_datExpr, AD_multiExpr[[i]]$datExpr)
AD_multiExpr[[i]]$datMeta$Subject = rownames(AD_multiExpr[[i]]$datMeta)
names(AD_multiExpr[[i]]$datMeta) <- gsub("disease_status","Group",names(AD_multiExpr[[i]]$datMeta))
AD_multiExpr[[i]]$datMeta$Group <- gsub("Control","CTL",AD_multiExpr[[i]]$datMeta$Group)
AD_multiExpr[[i]]$datMeta$Group <- gsub("Affected","AD",AD_multiExpr[[i]]$datMeta$Group)
AD_multiExpr[[i]]$datMeta$Sex <- AD_multiExpr[[i]]$datMeta$Sex
AD_datMeta = rbind(AD_datMeta, AD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])

i=2;
AD_datExpr = cbind(AD_datExpr, AD_multiExpr[[i]]$datExpr)
AD_multiExpr[[i]]$datMeta$Subject = rownames(AD_multiExpr[[i]]$datMeta)
names(AD_multiExpr[[i]]$datMeta) <- gsub("Disease_State","Group",names(AD_multiExpr[[i]]$datMeta))
AD_multiExpr[[i]]$datMeta$Group <- gsub("normal","CTL",AD_multiExpr[[i]]$datMeta$Group)
AD_multiExpr[[i]]$datMeta$Group <- gsub("Alzheimer's_Disease","AD",AD_multiExpr[[i]]$datMeta$Group)
AD_multiExpr[[i]]$datMeta$Sex <- AD_multiExpr[[i]]$datMeta$Sex

AD_datMeta = rbind(AD_datMeta, AD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])



AD_datMeta = AD_datMeta[-1,]
to_keep = AD_datMeta$Group %in% c("CTL", "AD")
AD_datExpr = AD_datExpr[,to_keep]; AD_datMeta = AD_datMeta[to_keep,]
AD_datMeta$Group = factor(AD_datMeta$Group, levels=c("CTL", "AD"))
AD_datMeta$Study = as.factor(AD_datMeta$Study)
AD_datMeta$Subject = as.factor(AD_datMeta$Subject)

AD_male_Meta <- AD_datMeta[AD_datMeta$Sex=="F"|AD_datMeta$Sex=="female",]
AD_male_expr <- AD_datExpr[,rownames(AD_male_Meta)]

male_AD_meta = matrix(NA, nrow=nrow(AD_male_expr), ncol=3)
for(i in 1:nrow(AD_male_expr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AD_male_expr[i,])
  tryCatch({
    male_AD_meta[i,] = summary(lme(expr~ Group + Study,data = AD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

male_AD_meta=as.data.frame(male_AD_meta)
colnames(male_AD_meta) = c("beta", "SE", "p")
rownames(male_AD_meta) = genes
male_AD_meta$fdr = p.adjust(male_AD_meta$p, "fdr")
male_AD_meta$symbol=AD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, AD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#AD_meta= AD_meta[!apply(is.na(AD_meta),1,any),]
write.csv(file="./results/tables/Microarray_AD_immunegene_female_metaanalysis_190422.csv",male_AD_meta)

#ASD
ASD_multiExpr = vector(mode="list",length = 3)
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned//Microarray_ASD_voineagu_normalized_CR_cleaned.rdata")
i=1; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_chow_normalized_CR_cleaned.rdata")
i=2; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_Garbett_normalized_CR_regressed.RData.rdata")
i=3; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)

genes= intersect(intersect(intersect(rownames(ASD_multiExpr[[1]]$datExpr), rownames(ASD_multiExpr[[2]]$datExpr)), rownames(ASD_multiExpr[[3]]$datExpr)),rownames(immunegene))

for(i in 1:3) ASD_multiExpr[[i]]$datExpr = ASD_multiExpr[[i]]$datExpr[match(genes,rownames(ASD_multiExpr[[i]]$datExpr)),]

ASD_datExpr = cbind(ASD_multiExpr[[1]]$datExpr,ASD_multiExpr[[2]]$datExpr,ASD_multiExpr[[3]]$datExpr)
ASD_datMeta = rbind(ASD_multiExpr[[1]]$datMeta[,c("Study", "Subject","Group","Sex")],ASD_multiExpr[[2]]$datMeta[,c("Study", "Subject","Group","Sex")],ASD_multiExpr[[3]]$datMeta[,c("Study", "Subject","Group","Sex")])

ASD_male_Meta <- ASD_datMeta[ASD_datMeta$Sex=="F"|ASD_datMeta$Sex=="female",]
ASD_male_expr <- ASD_datExpr[,rownames(ASD_male_Meta)]

asd_meta = matrix(NA, nrow=nrow(ASD_male_expr), ncol=3)
for(i in 1:nrow(ASD_male_expr)) {
  if(i%%100==0) print(i)
  expr = ASD_male_expr[i,]
  tryCatch({
    asd_meta[i,] = summary(lme(expr~ Group + Study,data = ASD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
asd_meta = as.data.frame(asd_meta)
colnames(asd_meta) = c("beta", "SE", "p")
rownames(asd_meta) = genes
asd_meta$fdr = p.adjust(asd_meta$p, "fdr")
asd_meta$symbol=ASD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, ASD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#to_keep = !apply(is.na(asd_meta),1,any)
#asd_meta = asd_meta[to_keep,]
write.csv(file="./results/tables//Microarray_ASD_immunegene_female_metaanalysis_190422.csv",asd_meta)

#SCZ
files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("SCZ",files)]
n=length(files)
SCZ_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  SCZ_multiExpr[[i]]$datExpr = datExpr 
  SCZ_multiExpr[[i]]$datMeta = datMeta
  SCZ_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(intersect(rownames(SCZ_multiExpr[[1]]$datExpr), rownames(SCZ_multiExpr[[2]]$datExpr)),rownames(immunegene))

for(i in 3:n) genes = intersect(genes, rownames(SCZ_multiExpr[[i]]$datExpr))
for(i in 1:n) SCZ_multiExpr[[i]]$datExpr = SCZ_multiExpr[[i]]$datExpr[match(genes,rownames(SCZ_multiExpr[[i]]$datExpr)),]

SCZ_datExpr = data.frame(row.names = genes)
SCZ_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)
for(i in 1:n) {
  SCZ_datExpr = cbind(SCZ_datExpr, SCZ_multiExpr[[i]]$datExpr)
  
  #if("Group.SCZ" %in% colnames(SCZ_multiExpr[[i]]$datMeta)) SCZ_multiExpr[[i]]$datMeta$Group = SCZ_multiExpr[[i]]$datMeta$Group.SCZ   #Only for unique controls
  SCZ_datMeta = rbind(SCZ_datMeta, SCZ_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])
  
}
SCZ_datMeta = SCZ_datMeta[-1,]
to_keep = SCZ_datMeta$Group %in% c("CTL", "SCZ")
SCZ_datExpr = SCZ_datExpr[,to_keep]; SCZ_datMeta = SCZ_datMeta[to_keep,]
SCZ_datMeta$Group = factor(SCZ_datMeta$Group)
SCZ_datMeta$Study = as.factor(SCZ_datMeta$Study)
SCZ_datMeta$Subject = as.factor(SCZ_datMeta$Subject)

SCZ_male_Meta <- SCZ_datMeta[SCZ_datMeta$Sex=="M"|SCZ_datMeta$Sex=="male",]
SCZ_male_expr <- SCZ_datExpr[,rownames(SCZ_male_Meta)]


scz_meta = matrix(NA, nrow=nrow(SCZ_datExpr), ncol=3)
for(i in 1:nrow(SCZ_male_expr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(SCZ_male_expr[i,])
  tryCatch({
    scz_meta[i,] = summary(lme(expr~ Group + Study,data = SCZ_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

scz_meta=as.data.frame(scz_meta)
colnames(scz_meta) = c("beta", "SE", "p")
rownames(scz_meta) = genes
scz_meta$fdr = p.adjust(scz_meta$p, "fdr")
scz_meta$symbol=SCZ_multiExpr[[1]]$datProbes$external_gene_id[match(genes, SCZ_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#scz_meta= scz_meta[!apply(is.na(scz_meta),1,any),]

write.csv(file="./results/tables/Microarray_SCZ_immunegene_male_metaanalysis_190422.csv",scz_meta)



#BD
files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("_BD_",files)]
n=length(files)
BD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  BD_multiExpr[[i]]$datExpr = datExpr 
  BD_multiExpr[[i]]$datMeta = datMeta
  BD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(intersect(rownames(BD_multiExpr[[1]]$datExpr), rownames(BD_multiExpr[[2]]$datExpr)),rownames(immunegene))
for(i in 3:n) genes = intersect(genes, rownames(BD_multiExpr[[i]]$datExpr))
for(i in 1:n) BD_multiExpr[[i]]$datExpr = BD_multiExpr[[i]]$datExpr[match(genes,rownames(BD_multiExpr[[i]]$datExpr)),]

BD_datExpr = data.frame(row.names = genes)
BD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)
for(i in 1:n) {
  BD_datExpr = cbind(BD_datExpr, BD_multiExpr[[i]]$datExpr)
  #if("Group.BD" %in% colnames(BD_multiExpr[[i]]$datMeta)) BD_multiExpr[[i]]$datMeta$Group = BD_multiExpr[[i]]$datMeta$Group.BD ## --> only for unique controls
  BD_datMeta = rbind(BD_datMeta, BD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])
}
BD_datMeta = BD_datMeta[-1,]
to_keep = BD_datMeta$Group %in% c("CTL", "BD")
BD_datExpr = BD_datExpr[,to_keep]; BD_datMeta = BD_datMeta[to_keep,]
BD_datMeta$Group = factor(BD_datMeta$Group, levels=c("CTL", "BD"))
BD_datMeta$Study = as.factor(BD_datMeta$Study)
BD_datMeta$Subject = as.factor(BD_datMeta$Subject)

BD_male_Meta <- BD_datMeta[BD_datMeta$Sex=="F"|BD_datMeta$Sex=="female",]
BD_male_expr <- BD_datExpr[,rownames(BD_male_Meta)]

bd_meta = matrix(NA, nrow=nrow(BD_datExpr), ncol=3)
for(i in 1:nrow(BD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(BD_datExpr[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ Group + Study,data = BD_datMeta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "p")
rownames(bd_meta) = genes
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
bd_meta$symbol=BD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, BD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#bd_meta= bd_meta[!apply(is.na(bd_meta),1,any),]

write.csv(file="./results/tables/Microarray_BD_immunegene_female_metaanalysis_190422.csv",bd_meta)


#MDD
files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("_MDD_",files)]
n=length(files)
MDD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  MDD_multiExpr[[i]]$datExpr = datExpr 
  MDD_multiExpr[[i]]$datMeta = datMeta
  MDD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(intersect(rownames(MDD_multiExpr[[1]]$datExpr), rownames(MDD_multiExpr[[2]]$datExpr)),rownames(immunegene))
for(i in 1:n) MDD_multiExpr[[i]]$datExpr = MDD_multiExpr[[i]]$datExpr[match(genes,rownames(MDD_multiExpr[[i]]$datExpr)),]

MDD_datExpr = data.frame(row.names = genes)
MDD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)
for(i in 1:n) {
  MDD_datExpr = cbind(MDD_datExpr, MDD_multiExpr[[i]]$datExpr)
  MDD_datMeta = rbind(MDD_datMeta, MDD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])
}
MDD_datMeta = MDD_datMeta[-1,]
to_keep = MDD_datMeta$Group %in% c("CTL", "MDD")
MDD_datExpr = MDD_datExpr[,to_keep]; MDD_datMeta = MDD_datMeta[to_keep,]
MDD_datMeta$Group = factor(MDD_datMeta$Group, levels=c("CTL", "MDD"))
MDD_datMeta$Study = as.factor(MDD_datMeta$Study)
MDD_datMeta$Subject = as.factor(MDD_datMeta$Subject)

MDD_male_Meta <- MDD_datMeta[MDD_datMeta$Sex=="F"|MDD_datMeta$Sex=="female",]
MDD_male_expr <- MDD_datExpr[,rownames(MDD_male_Meta)]

mdd_meta = matrix(NA, nrow=nrow(MDD_datExpr), ncol=3)
for(i in 1:nrow(MDD_male_expr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(MDD_male_expr[i,])
  tryCatch({
    mdd_meta[i,] = summary(lme(expr~ Group + Study,data = MDD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

mdd_meta=as.data.frame(mdd_meta)
colnames(mdd_meta) = c("beta", "SE", "p")
rownames(mdd_meta) = genes
mdd_meta$fdr = p.adjust(mdd_meta$p, "fdr")
mdd_meta$symbol=MDD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, MDD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#mdd_meta= mdd_meta[!apply(is.na(mdd_meta),1,any),]
write.csv(file="./results/tables/Microarray_MDD_immunegene_female_metaanalysis_190422.csv",mdd_meta)


##AAD
load("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AAD_mayfield_normalized_CR_regressed.RData")
AAD_datExpr = datExpr; rm(datExpr)
AAD_datMeta = datMeta; rm(datMeta)
AAD_datProbes = datProbes; rm(datProbes)

genes= intersect(rownames(AAD_datExpr),rownames(immunegene))
AAD_datExpr = AAD_datExpr[match(genes,rownames(AAD_datExpr)),]
AAD_Group = factor(AAD_datMeta$Group, levels=c("CTL", "ETOH"))
AAD_Study = as.factor(AAD_datMeta$Study)
AAD_Subject = as.factor(AAD_datMeta$Subject)
AAD_Sex = as.factor(AAD_datMeta$Sex)
AAD_datMeta=data.frame(Study=AAD_Study, Subject=AAD_Subject, Group=AAD_Group,Sex=AAD_Sex)


AAD_male_Meta <- AAD_datMeta[AAD_datMeta$Sex=="M"|AAD_datMeta$Sex=="male",]
AAD_male_expr <- AAD_datExpr[,rownames(AAD_male_Meta)]


aad_meta = matrix(NA, nrow=nrow(AAD_datExpr), ncol=3)
for(i in 1:nrow(AAD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AAD_datExpr[i,])
  tryCatch({
    aad_meta[i,] = summary(lme(expr~ Group,data = AAD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}
aad_meta=as.data.frame(aad_meta)
colnames(aad_meta) = c("beta", "SE", "p")
rownames(aad_meta) = genes
aad_meta$fdr = p.adjust(aad_meta$p, "fdr")
aad_meta$symbol=AAD_datProbes$`Associated Gene Name`[match(genes, AAD_datProbes$`Ensembl Gene ID`)]
#aad_meta= aad_meta[!apply(is.na(aad_meta),1,any),]
write.csv(file="./results/tables/Microarray_AAD_immuneGene_metaanalysis_123018.csv",aad_meta)


#PD
files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("PD",files)]
n=length(files)
PD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  PD_multiExpr[[i]]$datExpr = datExpr 
  PD_multiExpr[[i]]$datMeta = datMeta
  PD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(intersect(rownames(PD_multiExpr[[1]]$datExpr), rownames(PD_multiExpr[[2]]$datExpr)),rownames(immunegene))
for(i in 3:n) genes = intersect(genes, rownames(PD_multiExpr[[i]]$datExpr))
for(i in 1:n) PD_multiExpr[[i]]$datExpr = PD_multiExpr[[i]]$datExpr[match(genes,rownames(PD_multiExpr[[i]]$datExpr)),]
#names(PD_multiExpr[[1]]$datMeta) <- gsub("disease_state","Group",names(PD_multiExpr[[1]]$datMeta))
names(PD_multiExpr[[1]]$datMeta) <- gsub("disease","Group",names(PD_multiExpr[[1]]$datMeta))
names(PD_multiExpr[[2]]$datMeta) <- gsub("disease_state","Group",names(PD_multiExpr[[2]]$datMeta))
names(PD_multiExpr[[3]]$datMeta) <- gsub("DiseaseState","Group",names(PD_multiExpr[[3]]$datMeta))
PD_multiExpr[[2]]$datMeta$Source <- rownames(PD_multiExpr[[2]]$datMeta)

PD_datExpr = data.frame(row.names = genes)
PD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)
for(i in 1:n) {
  PD_datExpr = cbind(PD_datExpr, PD_multiExpr[[i]]$datExpr)
  
  #if("Group.PD" %in% colnames(PD_multiExpr[[i]]$datMeta)) PD_multiExpr[[i]]$datMeta$Group = PD_multiExpr[[i]]$datMeta$Group.PD #Only for unique controls
  names(PD_multiExpr[[i]]$datMeta) <- gsub("Source","Subject",names(PD_multiExpr[[i]]$datMeta))
  PD_datMeta = rbind(PD_datMeta, PD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])
  PD_datMeta$Group <- gsub("control","CTL",PD_datMeta$Group)
  PD_datMeta$Group <- gsub("normal","CTL",PD_datMeta$Group)
  PD_datMeta$Group <- gsub("Parkinson's_disease","PD",PD_datMeta$Group)
}
PD_datMeta = PD_datMeta[-1,]
to_keep = PD_datMeta$Group %in% c("CTL", "PD")
PD_datMeta$Group = factor(PD_datMeta$Group)
PD_datMeta$Study = as.factor(PD_datMeta$Study)
PD_datMeta$Subject = as.factor(PD_datMeta$Subject)

PD_male_Meta <- PD_datMeta[PD_datMeta$Sex=="F"|PD_datMeta$Sex=="female",]
PD_male_expr <- PD_datExpr[,colnames(PD_datExpr) %in% rownames(PD_male_Meta)]

PD_meta = matrix(NA, nrow=nrow(PD_male_expr), ncol=3)
for(i in 1:nrow(PD_male_expr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(PD_male_expr[i,])
  tryCatch({
    PD_meta[i,] = summary(lme(expr~ Group + Study,data = PD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

PD_meta=as.data.frame(PD_meta)
colnames(PD_meta) = c("beta", "SE", "p")
rownames(PD_meta) = genes
PD_meta$fdr = p.adjust(PD_meta$p, "fdr")
PD_meta$symbol=PD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, PD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#PD_meta= PD_meta[!apply(is.na(PD_meta),1,any),]

write.csv(file="./results/tables/Microarray_PD_immunegene_female_metaanalysis_190422.csv",PD_meta)


#IBD
files = dir("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("_IBD_",files)]
n=length(files)
IBD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("./working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  IBD_multiExpr[[i]]$datExpr = datExpr 
  IBD_multiExpr[[i]]$datMeta = datMeta
  IBD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(intersect(rownames(IBD_multiExpr[[1]]$datExpr), rownames(IBD_multiExpr[[2]]$datExpr)),rownames(immunegene))
for(i in 1:n) IBD_multiExpr[[i]]$datExpr = IBD_multiExpr[[i]]$datExpr[match(genes,rownames(IBD_multiExpr[[i]]$datExpr)),]

IBD_datExpr = data.frame(row.names = genes)
IBD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA,Sex=NA)
for(i in 1:n) {
  IBD_datExpr = cbind(IBD_datExpr, IBD_multiExpr[[i]]$datExpr)
  IBD_multiExpr[[i]]$datMeta$Subject = rownames(IBD_multiExpr[[i]]$datMeta)
  IBD_datMeta = rbind(IBD_datMeta, IBD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group","Sex")])
}
IBD_datMeta = IBD_datMeta[-1,]
to_keep = IBD_datMeta$Group %in% c("CTL", "IBD")
IBD_datExpr = IBD_datExpr[,to_keep]; IBD_datMeta = IBD_datMeta[to_keep,]
IBD_datMeta$Group = factor(IBD_datMeta$Group, levels=c("CTL", "IBD"))
IBD_datMeta$Study = as.factor(IBD_datMeta$Study)
IBD_datMeta$Subject = as.factor(IBD_datMeta$Subject)

IBD_male_Meta <- IBD_datMeta[IBD_datMeta$Sex=="M"|IBD_datMeta$Sex=="male",]
IBD_male_expr <- IBD_datExpr[,colnames(IBD_datExpr) %in% rownames(IBD_male_Meta)]

ibd_meta = matrix(NA, nrow=nrow(IBD_datExpr), ncol=3)
for(i in 1:nrow(IBD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(IBD_male_expr[i,])
  tryCatch({
    ibd_meta[i,] = summary(lme(expr~ Group + Study,data = IBD_male_Meta, random=~1|Subject))$tTable[2,c(1,2,5)]
  }, error=function(e){})
}

ibd_meta=as.data.frame(ibd_meta)
colnames(ibd_meta) = c("beta", "SE", "p")
rownames(ibd_meta) = genes
ibd_meta$fdr = p.adjust(ibd_meta$p, "fdr")
ibd_meta$symbol=IBD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, IBD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
#ibd_meta= ibd_meta[!apply(is.na(ibd_meta),1,any),]
write.csv(file="./results/tables/Microarray_IBD_immunegene_male_metaanalysis_123018.csv",ibd_meta)


##Save compiled expression and metadata for permutation testing
multiExpr = vector(mode="list", length=8)
multiExpr[[1]]$datExpr = ASD_datExpr;  multiExpr[[1]]$datMeta = ASD_datMeta 
multiExpr[[2]]$datExpr = SCZ_datExpr; multiExpr[[2]]$datMeta = SCZ_datMeta 
multiExpr[[3]]$datExpr = BD_datExpr; multiExpr[[3]]$datMeta = BD_datMeta 
multiExpr[[4]]$datExpr = MDD_datExpr; multiExpr[[4]]$datMeta = MDD_datMeta 
multiExpr[[5]]$datExpr = AAD_datExpr; multiExpr[[5]]$datMeta = AAD_datMeta 
multiExpr[[6]]$datExpr = IBD_datExpr; multiExpr[[6]]$datMeta = IBD_datMeta
multiExpr[[7]]$datExpr = AD_datExpr; multiExpr[[7]]$datMeta = AD_datMeta
multiExpr[[8]]$datExpr = PD_datExpr; multiExpr[[8]]$datMeta = PD_datMeta
names(multiExpr) = c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD","AD","PD")

genes = intersect(intersect(rownames(multiExpr[[1]]$datExpr), rownames(multiExpr[[2]]$datExpr)),rownames(immunegene))
for(i in 3:8) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))

for(i in 1:8) multiExpr[[i]]$datExpr = multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),]
save(file="./working_data/Microarray/Microarray_immunegene_compiledForPermutationTesting123018.RData", multiExpr)

