# Generate forest plots with the results of the logistic regression analysis per cancer type.
# library(tidyverse)

# Load results 
logreg_per_cancer_wp<- readRDS("results/data/logreg_wp.rds")
logreg_per_cancer_wp$cancer_type<- rownames(logreg_per_cancer_wp)
logreg_per_cancer_wm<- readRDS("results/data/logreg_wm.rds")
logreg_per_cancer_wm$cancer_type<- rownames(logreg_per_cancer_wm)

# Pan cancer analysis
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_mut_HLA_panCancer_analysis.pdf")
  else svglite::svglite("results/figs/manuscript_mut_HLA_panCancer_analysis.svg")
  par(mfrow=c(1,2))
  bp<- barplot(c(logreg_per_cancer_wm["PAN","OR"],logreg_per_cancer_wp["PAN","OR"]),ylim=c(0.95,1.2),xpd=F,ylab="Odds ratio",names.arg = c("Within mutations","Within patients"))
  arrows(bp[,1],c(logreg_per_cancer_wm["PAN","ci_low"],logreg_per_cancer_wp["PAN","ci_low"]),bp[,1],c(logreg_per_cancer_wm["PAN","ci_high"],logreg_per_cancer_wp["PAN","ci_high"]),angle = 90, code = 3,length = 0.1,xpd=NA)
  abline(h=1, lty=2)
  text(bp, rep(1.22,2),labels = signif(c(logreg_per_cancer_wm["PAN","p"],logreg_per_cancer_wp["PAN","p"]),3),xpd=NA)
  dev.off()
}

# Per cancer analysis
logreg_per_cancer_wp<- logreg_per_cancer_wp[logreg_per_cancer_wp$cancer_type!="PAN",]
logreg_per_cancer_wp$q<- p.adjust(logreg_per_cancer_wp$p,"fdr")
logreg_per_cancer_wm<- logreg_per_cancer_wm[logreg_per_cancer_wm$cancer_type!="PAN",]
logreg_per_cancer_wm$q<- p.adjust(logreg_per_cancer_wm$p,"fdr")

# Plot per cancer
source("scripts/functions/plot_OR.R")
library(gridExtra)
wm<- plot_OR(logreg_per_cancer_wm,"Within mutations")
wp<- plot_OR(logreg_per_cancer_wp,"Within patients")

for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_mut_HLA_perCancer_analysis.pdf")
  else svglite::svglite("results/figs/manuscript_mut_HLA_perCancer_analysis.svg",height = 10)
  grid.arrange(a=wm,wp,ncol=2)
  dev.off()
}

