# Load data
###########

# Load mut matrix
mut_matrix<- readRDS(file="data/mut_matrix.rds")

# Load VP HLA matrices
HLA_matrix<- readRDS(file="data/HLA_matrix.rds")

for(i in 1:2){
  if(i==1) pdf("results/figs/TCGA_analysis.pdf")

  # 1 PHBR ~ Frequency
  ######################
  mut_freq<- 100*colMeans(mut_matrix)
  # Use median PHBR per mutation
  PHBR_med<- apply(HLA_matrix,2,"median",na.rm=T)
  PHBR_freq_t<- as.data.frame(na.omit(cbind(mut_freq,PHBR=PHBR_med[names(mut_freq)])))
  
  # Plot
  if(i==2) svglite::svglite("results/figs/manuscript_TCGA_analysis_mutFreq.svg")
  plot(PHBR_freq_t$mut_freq,PHBR_freq_t$PHBR,pch=16,col="blue",frame.plot = F,ylim=c(0,12),xlab="Mutation frequency (%)",ylab="PHBR",xpd=NA)
  points(loess.smooth(PHBR_freq_t$mut_freq,PHBR_freq_t$PHBR,span = 1),type="l",lwd=3,col="red")
  abline(h=median(PHBR_freq_t$PHBR),lty=2)
  
  # Show 13 mutations
  muts<- readRDS("data/driver_muts13.rds")
  text(PHBR_freq_t[muts,"mut_freq"],PHBR_freq_t[muts,"PHBR"],muts,cex=0.7,col="red")
  points(PHBR_freq_t[muts,"mut_freq"],PHBR_freq_t[muts,"PHBR"],pch=16,col="red")
  
  # Show other frequent mutations
  mut_freq<- setdiff(rownames(PHBR_freq_t[PHBR_freq_t$mut_freq>1,]),muts)
  text(PHBR_freq_t[mut_freq,"mut_freq"],PHBR_freq_t[mut_freq,"PHBR"],mut_freq,cex=0.7,col="black")
  points(PHBR_freq_t[mut_freq,"mut_freq"],PHBR_freq_t[mut_freq,"PHBR"],pch=16,col="black")
  
  if(i==2) dev.off()
  
  # 2 Focus on top 10
  ####################
  if(i==2) svglite::svglite("results/figs/manuscript_TCGA_analysis_boxplot_top10.svg")
  par(mfrow=c(1,2))
  boxplot(list(PHBR_med[1:10],PHBR_med[1:30],PHBR_med[1:100],PHBR_med),outline=F,frame.plot=F,names=c("Top 10", "Top 30", "Top 100","All"),las=2,ylab="PHBR")
  if(i==2) dev.off()

  # Stats?
  summary(PHBR_med[1:10])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.950   2.830   3.727   3.855   4.523   6.181 
  summary(PHBR_med)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #  0.0714  1.0390  1.9393  2.8506  3.4664 20.4044       3 
  
  # Cumulative frequency?
  sum(colSums(mut_matrix,na.rm=T)[1:10])/sum(colSums(mut_matrix,na.rm=T)) #  0.2550916
  sum(colSums(mut_matrix,na.rm=T)[1:30])/sum(colSums(mut_matrix,na.rm=T)) # 0.3898167
  sum(colSums(mut_matrix,na.rm=T)[1:100])/sum(colSums(mut_matrix,na.rm=T)) #  0.5643585
  
  # 3 LOO
  ########
  
  # See scripts/manuscript_LOO_analysis_TCGA.R
  
  # Baseline OR
  baseline<- as.data.frame(readRDS("results/data/TCGA_logreg_wp.rds"))["mut",]
  
  # # OR after exclusion 13 muts
  # OR_excl13<- readRDS("results/data/vp_logreg_loo_mut13.rds")
  # 
  # LOO
  loo<- readRDS("results/data/vp_logreg_loo_TCGA.rds")
  loo<- loo[order(loo$OR),]
  
  # Calculate % contribution based on LOO
  loo$rel_effect<- 1-(as.numeric(loo$OR) - 1)/(baseline$OR-1)
  
  # Get muts
  muts<- readRDS("data/driver_muts13.rds")
  
  # Barplot
  if(i==2) svglite::svglite("results/figs/manuscript_TCGA_analysis_loo.svg")
  barplot(100*loo$rel_effect[1:20],names.arg = loo$mut[1:20],las=2,ylab= "Effect (%)",col=loo$mut%in%muts)
  abline(h=1,lty=2)
  if(i==2)   dev.off()

  # 4 dPHBR (ECDF)
  ################
  
  HLA_wt_matrix<- readRDS(file="data/HLA_wt_matrix.rds")
  dPHBR_matrix<- HLA_matrix-HLA_wt_matrix
  dPHBR_df<- data.frame(dPHBR=apply(dPHBR_matrix,2,"median",na.rm=T),row.names = colnames(dPHBR_matrix))
  
  if(i==2) svglite::svglite("results/figs/manuscript_TCGA_analysis_dPHBR_ecdf.svg")
  par(mfrow=c(1,2))
  dPHBR_ecdf<- ecdf(dPHBR_df$dPHBR)
  plot(dPHBR_ecdf,frame.plot=F,xlab="dPHBR",ylab = "Cumulative frequency")
  for(mut in muts){
    mut_co<- c(dPHBR_df[mut,"dPHBR"],dPHBR_ecdf(dPHBR_df[mut,"dPHBR"]))
    points(mut_co[1],mut_co[2],pch=16,col="red")
    text(mut_co[1],mut_co[2],gsub("_"," ",mut),cex=0.5,ad=c(0,1),col="red")
  }
  if(i==2) dev.off()
  cat(muts[1],c(dPHBR_df[muts[1],"dPHBR"],dPHBR_ecdf(dPHBR_df[muts[1],"dPHBR"])))
  cat(muts[11],c(dPHBR_df[muts[11],"dPHBR"],dPHBR_ecdf(dPHBR_df[muts[11],"dPHBR"])))
  # BRAF_V600E 2.703973 0.9780059
  # TP53_Y220C 3.548239 0.9868035
  
  # 5 PHBR_wt
  #############
  
  # OR
  logreg_wp_mut<- as.data.frame(readRDS("results/data/TCGA_logreg_wp.rds"))
  logreg_wp_mut$OR 
  # [1] 1.193 1.1056
  (as.numeric(logreg_wp_mut["wt","OR"])-1)/(as.numeric(logreg_wp_mut["mut","OR"])-1) # 54.7%
  # [1] 0.54688
  
  # Plot
  if(i==2) svglite::svglite("results/figs/manuscript_TCGA_analysis_OR_wt.svg")
  par(mfrow=c(1,2))
  bp<- barplot(logreg_wp_mut$OR,ylim=c(0.95,1.21),xpd=F,ylab="Odds ratio",names.arg = c("Mut","WT"))
  arrows(bp[,1],logreg_wp_mut$`2.5 %`,bp[,1],logreg_wp_mut$`97.5 %`,angle = 90, code = 3,length = 0.1,xpd=NA)
  abline(h=1, lty=2)
  text(bp, rep(1.22,2),labels = signif(logreg_wp_mut$p,3),xpd=NA)
  dev.off()
}  

