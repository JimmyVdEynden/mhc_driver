# Plot
######
pop_PHBR_ls<- readRDS("results/data/1kg_LR.rds")
pop_PHBR_ls_lowFreq<- readRDS("results/data/1kg_LR_lowFreq.rds")
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")
muts<- readRDS("data/driver_muts13.rds")

# Boxplots
##########
# pdf("results/figs/1kg_analysis_boxplot.pdf",height = 7, width = 7)
svglite::svglite("results/figs/1kg_analysis_boxplot.svg",height = 7, width = 7)
par(mfrow=c(2,5))

for(i in 1:2){
  
  if(i==1) pop_PHBR<- pop_PHBR_ls
  else pop_PHBR<- pop_PHBR_ls_lowFreq
  
  for(pop in names(pop_PHBR_ls)){
    
    vp_HLA1_mut_matrix<- pop_PHBR[[pop]][["vp_HLA1_mut_matrix"]]
    
    # Vectors
    vp_HLA_vector<- as.numeric(vp_HLA1_mut_matrix)
    vp_mut_vector<- as.logical(vp_mut_matrix!=0)
    vp_mut_vector<- vp_mut_vector[!is.na(vp_HLA_vector)]
    vp_HLA_vector<- vp_HLA_vector[!is.na(vp_HLA_vector)]
    
    # Mutations lower PHBR score? 
    cat("Obs: ",median(vp_HLA_vector[vp_mut_vector==T],na.rm = T),"\n") 
    cat("Unobs: ",median(vp_HLA_vector[vp_mut_vector==F],na.rm = T),"\n")
    p_tmp<-  wilcox.test(vp_HLA_vector~vp_mut_vector)$p.value 
    
    bp<- boxplot(vp_HLA_vector~vp_mut_vector,outline=F,names=c("-","+"),xlab="Mutated",ylab="PHBR",frame.plot=F,main=pop)
    text(1.5,max(bp$stats),labels = paste0("P=",format(p_tmp,scientific = T,digits=3)),xpd=NA)
  }
}
dev.off()

# Logreg
#########
# pdf("results/figs/1kg_analysis_LR.pdf",height = 7, width = 7)
svglite::svglite("results/figs/1kg_analysis_LR.svg",height = 7, width = 7)
par(mfrow=c(2,2))

logreg_LR<- t(sapply(pop_PHBR_ls, function(x) x[["logreg"]]["wp",]))                             
bp<- barplot(logreg_LR[,"OR"],ylim=c(0.95,1.25),ylab="Odds ratio",xpd=F)
arrows(bp[,1],logreg_LR[,"2.5 %"],bp[,1],logreg_LR[,"97.5 %"],angle = 90, code = 3,length = 0.1,xpd=NA)
abline(h=1, lty=2)
text(bp, rep(1.27,2),labels = signif(logreg_LR[,"p"],3),xpd=NA)

logreg_LR<- t(sapply(pop_PHBR_ls_lowFreq, function(x) x[["logreg"]]["wp",]))                             
bp<- barplot(logreg_LR[,"OR"],ylim=c(0.95,1.25),ylab="Odds ratio",xpd=F)
arrows(bp[,1],logreg_LR[,"2.5 %"],bp[,1],logreg_LR[,"97.5 %"],angle = 90, code = 3,length = 0.1,xpd=NA)
abline(h=1, lty=2)
text(bp, rep(1.27,2),labels = signif(logreg_LR[,"p"],3),xpd=NA)

dev.off()

