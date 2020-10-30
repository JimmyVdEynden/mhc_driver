############################
# Plot LOO results
############################

# Baseline OR
baseline<- as.data.frame(readRDS("results/data/vp_logreg_wp.rds"))["mut",]

# OR after exclusion 13 muts
OR_excl13<- readRDS("results/data/vp_logreg_loo_mut13.rds")

# LOO
loo<- readRDS("results/data/vp_logreg_loo.rds")
loo<- loo[order(loo$OR),]

# Calculate % contribution based on LOO
loo$rel_effect<- 1-(loo$OR - 1)/(baseline$OR-1)

# Get muts
muts<- readRDS("data/driver_muts13.rds")

# Barplot
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_vp_loo.pdf")
  else svglite::svglite("results/figs/manuscript_vp_loo.svg")
  barplot(100*loo$rel_effect[1:20],names.arg = loo$mutation[1:20],las=2,ylab= "Effect (%)")
  abline(h=1,lty=2)
  dev.off()
}

# Illustrate effect on BRAF + cumul
#####################################

BRAF_t<- data.frame(OR=c(baseline$OR,loo$OR[1],OR_excl13["OR"]),CI_low=c(baseline$`2.5 %`,loo$`2.5 %`[1],OR_excl13["2.5 %"]),CI_high=c(baseline$`97.5 %`,loo$`97.5 %`[1],OR_excl13["97.5 %"]))

for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_vp_loo_concept.pdf")
  else svglite::svglite("results/figs/manuscript_vp_loo_concept.svg")
  plot(BRAF_t$OR,c(1,2,3),ylim=c(0,4),frame.plot=F,axes=F,xlab="OR",ylab=NA,xlim=c(0.95,1.25))
  axis(1)
  for(i in 1:3) lines(c(BRAF_t$CI_low[i],BRAF_t$CI_high[i]),c(i,i))
  lines(c(1,baseline$OR),c(4,4))
  abline(v=c(1,loo$OR[1],baseline$OR),lty=2)
  dev.off()
}

