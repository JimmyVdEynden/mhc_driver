# Load data
###########

# Load VP mut matrix
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")

# Load VP HLA matrices
load(file="data/vp_HLA_matrix.RData")

# Analysis
##########
for(i in 1:2){
  if(i==1) pdf("results/figs/vp_boxplot.pdf")
  else svglite::svglite("results/figs/vp_boxplot.svg")
  par(mfrow=c(1,3))
  
  for(i in 1:3){
    
    cat("\n",i,"\n")
    vp_HLA_matrix<- list(vp_HLA1_mut_matrix,vp_HLA1_wt_matrix,vp_HLA2_mut_matrix)[[i]]
    
    # Get score vector
    vp_HLA_vector<- as.numeric(vp_HLA_matrix)
    vp_mut_vector<- as.logical(vp_mut_matrix!=0)
    vp_mut_vector<- vp_mut_vector[!is.na(vp_HLA_vector)]
    vp_HLA_vector<- vp_HLA_vector[!is.na(vp_HLA_vector)]
    
    # Mutations lower PHBR score? (Fig. 4b)
    cat("Obs: ",median(vp_HLA_vector[vp_mut_vector==T],na.rm = T),"\n") 
    cat("Unobs: ",median(vp_HLA_vector[vp_mut_vector==F],na.rm = T),"\n")
    p_tmp<-  wilcox.test(vp_HLA_vector~vp_mut_vector)$p.value 
    
    bp<- boxplot(vp_HLA_vector~vp_mut_vector,outline=F,names=c("-","+"),xlab="Mutated",ylab="PHBR",frame.plot=F)
    text(1.5,max(bp$stats),labels = paste0("P=",format(p_tmp,scientific = T,digits=3)),xpd=NA)
  }
  dev.off()
}

# 1 
# Obs:  2.023327 
# Unobs:  1.544127 
# 
# 2 
# Obs:  1.91078 
# Unobs:  1.506813 
#
# 3 
# Obs:  20.35998 
# Unobs:  17.33244 

# Show frequency effect
########################

# 1 PHBR ~ Frequency
######################
load("data/PHBR_proto.RData")

mut_freq<- 100*colSums(vp_mut_matrix)/nrow(vp_mut_matrix)
PHBR_freq_t<- as.data.frame(na.omit(cbind(mut_freq,PHBR=PHBR_proto_mhc1_mut[names(mut_freq)])))

# Plot
for(i in 1:2){
  if(i==1) pdf("results/figs/vp_PHBR_mutFreq.pdf")
  else svglite::svglite("results/figs/vp_PHBR_mutFreq.svg")
  plot(PHBR_freq_t$mut_freq,PHBR_freq_t$PHBR,pch=16,col="blue",frame.plot = F,ylim=c(0,12),xlab="Mutation frequency (%)",ylab="PHBR")
  points(loess.smooth(PHBR_freq_t$mut_freq,PHBR_freq_t$PHBR,span = 1),type="l",lwd=3,col="red")
  abline(h=median(PHBR_freq_t$PHBR),lty=2)
  
  # Show 13 mutations
  muts<- readRDS("data/driver_muts13.rds")
  text(PHBR_freq_t[muts,"mut_freq"],PHBR_freq_t[muts,"PHBR"],muts,cex=0.7)
  
  # Show loess after excl
  PHBR_freq_t_excl<- PHBR_freq_t[!rownames(PHBR_freq_t)%in%muts,]
  points(loess.smooth(PHBR_freq_t_excl$mut_freq,PHBR_freq_t_excl$PHBR,span = 1),type="l",lwd=3,col="red",lty=2)
  dev.off()
  }  

# 2 Focus on top 10
####################
for(i in 1:2){
  if(i==1) pdf("results/figs/vp_boxplot_top10.pdf")
  else svglite::svglite("results/figs/vp_boxplot_top10.svg")
  par(mfrow=c(1,2))
  boxplot(list(colMeans(vp_HLA1_mut_matrix,na.rm=T)[1:10],colMeans(vp_HLA1_mut_matrix,na.rm=T)[1:30],colMeans(vp_HLA1_mut_matrix,na.rm=T)[1:100],colMeans(vp_HLA1_mut_matrix,na.rm=T)),outline=F,frame.plot=F,names=c("Top 10", "Top 30", "Top 100","All"),las=2,ylab="PHBR")
  boxplot(list(colMeans(vp_HLA2_mut_matrix,na.rm=T)[1:10],colMeans(vp_HLA2_mut_matrix,na.rm=T)[1:30],colMeans(vp_HLA2_mut_matrix,na.rm=T)[1:100],colMeans(vp_HLA2_mut_matrix,na.rm=T)),outline=F,frame.plot=F,names=c("Top 10", "Top 30", "Top 100","All"),las=2,ylab="PHBR")
  dev.off()
}


# Stats?
summary(colMeans(vp_HLA1_mut_matrix,na.rm=T)[1:10])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.423   2.057   3.925   4.680   5.454  10.984 
summary(colMeans(vp_HLA1_mut_matrix,na.rm=T))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01931  0.56140  1.54413  2.85121  3.62629 25.11809        3 
summary(colMeans(vp_HLA2_mut_matrix,na.rm=T)[1:10])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.485  20.064  31.640  31.008  39.911  58.597 
summary(colMeans(vp_HLA2_mut_matrix,na.rm=T))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#  0.02846  7.42330 17.35895 22.22450 31.99583 89.85449        6 
 
# Cumulative frequency?
sum(colSums(vp_mut_matrix,na.rm=T)[1:10])/sum(colSums(vp_mut_matrix,na.rm=T)) #  0.2506391
sum(colSums(vp_mut_matrix,na.rm=T)[1:30])/sum(colSums(vp_mut_matrix,na.rm=T)) # 0.393803
sum(colSums(vp_mut_matrix,na.rm=T)[1:100])/sum(colSums(vp_mut_matrix,na.rm=T)) #  0.5672359

