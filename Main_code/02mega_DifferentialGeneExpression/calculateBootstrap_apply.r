##2b_calculateNullDistribution

seed = 20191117
set.seed(seed)

options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R"); biocLite("lme4")
library(lme4)
setwd("../../ImmuneGeneAnalysis/")
load("./working_data/Microarray/Microarray_compiledForPermutationTesting20190407.RData") ##Compiled expression & metadata from 2a
commandArgs(trailingOnly=TRUE)
genes = rownames(multiExpr[[1]]$datExpr)
allmeta = matrix(NA,nrow=length(genes), length(multiExpr))
colnames(allmeta) = c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD","AD","PD")
allmeta=  as.data.frame(allmeta)

c <- lmerControl(check.nobs.vs.nlev = "ignore",
            check.nobs.vs.rankZ = "ignore",
            check.nobs.vs.nRE="ignore")

lmer_apply=function(x, datMeta) {
  if(length(unique(datMeta$Study)) == 1){
    return(summary(lmer(x ~ Group  + (1 | Subject),data=datMeta,control = c))$coefficients[2,1])
  }
  else if(length(unique(datMeta$Subject)) < nrow(datMeta)) {
    return(summary(lmer(x ~ Group + Study + (1 | Subject),data=datMeta))$coefficients[2,1])
  } else if(length(unique(datMeta$Subject)) > 1) {
    return(summary(lm(x ~ Group + Study,data=datMeta))$coefficients[2,1])
  } else {
    return(summary(lm(x ~ Group + Study,data=datMeta))$coefficients[2,1])
  }
}

#lmer_apply=function(x, datMeta) {
  if(length(unique(datMeta$Subject)) < nrow(datMeta)) {
    return(summary(lmer(x ~ Group + Study + (1 | Subject),data=datMeta))$coefficients[2,1])
  } else if(length(unique(datMeta$Subject)) > 1) {
    return(summary(lm(x ~ Group + Study,data=datMeta))$coefficients[2,1])
  } else {
    return(summary(lm(x ~ Group + Study,data=datMeta))$coefficients[2,1])
  }
}


for(i in 1:length(multiExpr)) {
  print(i)
  tt = matrix(NA, nrow=length(genes), ncol=3)
  subj = unique(as.character(multiExpr[[i]]$datMeta$Subject))
  subj_group = data.frame(Subject = subj, Group = multiExpr[[i]]$datMeta$Group[match(subj, multiExpr[[i]]$datMeta$Subject)])
  subj_group$Group = subj_group$Group[order(runif(nrow(subj_group)))] ##Randomly shuffle group assignment for each subject
  multiExpr[[i]]$datMeta$Group = subj_group$Group[match(multiExpr[[i]]$datMeta$Subject,subj_group$Subject)]
  allmeta[,i] = apply(multiExpr[[i]]$datExpr,1,lmer_apply,multiExpr[[i]]$datMeta)
  #for (n in 1:length(genes)){
    #allmeta[i,n] = lmer_apply(multiExpr[[i]]$datExpr[n,],multiExpr[[i]]$datMeta)
  #}
}

cor_vec = vector(mode="numeric")
comparisons = t(combn(seq(1,ncol(allmeta)),2))

for(i in 1:nrow(comparisons)) {
  r = cor(allmeta[,comparisons[i,1]], allmeta[,comparisons[i,2]], method = "spearman", use="pairwise.complete.obs")
  cor_vec = c(cor_vec,r)
}

write.table(cor_vec, file=paste("./working_data/NullDistribution/",seed, ".txt",sep=""),row.names = F,col.names =F)
