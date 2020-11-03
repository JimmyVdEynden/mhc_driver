# Logreg of WT vs MUT peptides
###############################

# OR
vp_logreg_wp_mut<- as.data.frame(readRDS("results/data/vp_logreg_wp.rds"))

vp_logreg_wp_mut$OR 
# [1] 1.200097 1.099298

(as.numeric(vp_logreg_wp_mut["wt","OR"])-1)/(as.numeric(vp_logreg_wp_mut["mut","OR"])-1) # 49.6%
# [1] 0.4962491

# Plot
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_vp_OR_wt.pdf")
  else svglite::svglite("results/figs/manuscript_vp_OR_wt.svg")
  par(mfrow=c(1,2))
  bp<- barplot(vp_logreg_wp_mut$OR,ylim=c(0.95,1.21),xpd=F,ylab="Odds ratio",names.arg = c("Mut","WT"))
  arrows(bp[,1],vp_logreg_wp_mut$`2.5 %`,bp[,1],vp_logreg_wp_mut$`97.5 %`,angle = 90, code = 3,length = 0.1,xpd=NA)
  abline(h=1, lty=2)
  text(bp, rep(1.22,2),labels = signif(vp_logreg_wp_mut$p,3),xpd=NA)
  dev.off()
}


