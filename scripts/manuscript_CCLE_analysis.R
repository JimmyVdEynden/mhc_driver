# Load CCLE data
###################
load("data/CCLE_data.RData")
sample_info<- as.data.frame(readr::read_csv("downloads/CCLE/sample_info.csv"))

# Remove NA rows
idx_NA<- which(is.na(HLA_matrix[,1]))
HLA_matrix<- HLA_matrix[-idx_NA,]
mut_matrix<- mut_matrix[-idx_NA,] # 781 x 688

# Only sample info for analyzed mutations
sample_info<- sample_info[sample_info$DepMap_ID%in%rownames(mut_matrix),]

##########################
# Real genotype analysis
##########################

for(i in 1:2){
  if(i==1) pdf("results/figs/CCLE_analysis.pdf")
  
  # 1 Overview data!
  ###################
  CCLE_cancer_t<- table(sample_info$primary_disease)
  sample_info$primary_disease[sample_info$primary_disease%in%names(CCLE_cancer_t[CCLE_cancer_t<20])]<- "Other"
  CCLE_cancer_t<- table(sample_info$primary_disease)
  CCLE_cancer_t<- sort(CCLE_cancer_t)
  CCLE_cancer_t
  
  # Bladder Cancer              Kidney Cancer               Liver Cancer             Gastric Cancer 
  # 22                         22                         22                         23 
  # Myeloma Endometrial/Uterine Cancer       Head and Neck Cancer          Pancreatic Cancer 
  # 24                         25                         26                         33 
  # Fibroblast    Colon/Colorectal Cancer                Skin Cancer             Ovarian Cancer 
  # 35                         39                         39                         40 
  # Breast Cancer               Brain Cancer                   Lymphoma                   Leukemia 
  # 41                         43                         49                         71 
  # Other                Lung Cancer 
  # 82                        145 
  
  # 25-color palette
  c25<-c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black","gold1",
         "skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise",
         "green1","yellow4","yellow3","darkorange4","brown")
  names(c25)<- c('Colon/Colorectal Cancer','Bladder Cancer','Breast Cancer','CESC','Myeloma','COAD','Lymphoma','GBM','Head and Neck Cancer','TGCT','PCPG','Kidney Cancer','Leukemia','Brain Cancer','Liver Cancer','Lung Cancer','Ovarian Cancer','Pancreatic Cancer','PRAD','Fibroblast','Skin Cancer','Gastric Cancer','THYM','Endometrial/Uterine Cancer')
  
  if(i==2) svglite::svglite("results/figs/CCLE_overview.svg")
  pie(CCLE_cancer_t,col=c25[names(CCLE_cancer_t)])
  if(i==2) dev.off()
  
  # 2 Box plots
  ################

  # Vectors
  HLA_vector<- as.numeric(HLA_matrix)
  mut_vector<- as.logical(mut_matrix!=0)
  mut_vector<- mut_vector[!is.na(HLA_vector)]
  HLA_vector<- HLA_vector[!is.na(HLA_vector)]
  
  # Mutations lower PHBR score? 
  cat("Obs: ",median(HLA_vector[mut_vector==T],na.rm = T),"\n") # 2.566098
  cat("Unobs: ",median(HLA_vector[mut_vector==F],na.rm = T),"\n") # 1.736275 
  p_tmp<-  wilcox.test(HLA_vector~mut_vector)$p.value # 1.50 e-16
  
  if(i==2) svglite::svglite("results/figs/CCLE_boxplot.svg")
  par(mfrow=c(2,4))
  bp<- boxplot(HLA_vector~mut_vector,outline=F,names=c("-","+"),xlab="Mutated",ylab="PHBR",frame.plot=F)
  text(1.5,max(bp$stats),labels = paste0("P=",format(p_tmp,scientific = T,digits=3)),xpd=NA)
  if(i==2) dev.off()
  
  # 3 Logreg
  ############
  source("scripts/functions/do_logreg.R")
  logreg_wp<- do_logreg(mut_matrix,HLA_matrix,"wp") # 1.24, p=e-13
  logreg_wm<- do_logreg(mut_matrix,HLA_matrix,"wm") # 1.05, p= 0,17

  if(i==2) svglite::svglite("results/figs/CCLE_LR.svg")
  par(mfrow=c(2,4))
  bp<- barplot(c(logreg_wm["OR"],logreg_wp["OR"]),ylim=c(0.95,1.30),xpd=F,ylab="Odds ratio",names.arg = c("Within mutations","Within patients"))
  arrows(bp[,1],c(logreg_wm["2.5 %"],logreg_wp["2.5 %"]),bp[,1],c(logreg_wm["97.5 %"],logreg_wp["97.5 %"]),angle = 90, code = 3,length = 0.1,xpd=NA)
  abline(h=1, lty=2)
  text(bp, rep(1.27,2),labels = signif(c(logreg_wm["p"],logreg_wp["p"]),3),xpd=NA)
  if(i==2) dev.off()

  # 4 PHBR ~ Frequency
  ######################
  mut_freq<- 100*colMeans(mut_matrix)
  
  # Use median PHBR per mutation
  PHBR_med<- apply(HLA_matrix,2,"median",na.rm=T)
  PHBR_freq_t<- as.data.frame(na.omit(cbind(mut_freq,PHBR=PHBR_med[names(mut_freq)])))
  
  # Plot
  if(i==2) svglite::svglite("results/figs/CCLE_mutFreq.svg")
  par(mfrow=c(1,1))
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
  
  # 5 Focus on top 10
  ####################
  if(i==2) svglite::svglite("results/figs/CCLE_boxplot_top10.svg")
  par(mfrow=c(1,2))
  boxplot(list(PHBR_med[1:10],PHBR_med[1:30],PHBR_med[1:100],PHBR_med),outline=F,frame.plot=F,names=c("Top 10", "Top 30", "Top 100","All"),las=2,ylab="PHBR")
  dev.off()
  
  # Stats?
  summary(PHBR_med[1:10])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.053   2.588   3.361   3.791   4.342   7.923 
  summary(PHBR_med)
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
  # 0.03373  0.88353  1.73779  2.74982  3.38447 22.34491        3 
  
  # Cumulative frequency?
  sum(colSums(mut_matrix,na.rm=T)[1:10])/sum(colSums(mut_matrix,na.rm=T)) #  0.3270525
  sum(colSums(mut_matrix,na.rm=T)[1:30])/sum(colSums(mut_matrix,na.rm=T)) # 0.4899058
  sum(colSums(mut_matrix,na.rm=T)[1:100])/sum(colSums(mut_matrix,na.rm=T)) #  0.7362046
 }  


###########################
# Analysis prototypes
##########################

Pmut<- colMeans(mut_matrix)

# Create virtual pts
vp_mut_matrix<- matrix(NA,nrow(mut_matrix),ncol(mut_matrix))
colnames(vp_mut_matrix)<- colnames(mut_matrix)
for(i in 1:nrow(vp_mut_matrix)){
  cat(i," ")
  for(j in 1:ncol(vp_mut_matrix)){
    vp_mut_matrix[i,j]<- sample(x = c(0,1),size = 1,replace = T, prob = c(1-Pmut[j],Pmut[j]))
  }
}

load("data/PHBR_proto.RData")
# 1) Mut mhc1 peptides
vp_HLA1_mut_matrix<- vp_mut_matrix
vp_HLA1_mut_matrix[,]<- NA
common_muts<- intersect(names(PHBR_proto_mhc1_mut),colnames(vp_HLA1_mut_matrix))
for(i in 1:nrow(vp_HLA1_mut_matrix)){
  vp_HLA1_mut_matrix[i,common_muts]<- PHBR_proto_mhc1_mut[common_muts]
}

mut_matrix<- vp_mut_matrix
HLA_matrix<- vp_HLA1_mut_matrix
 
for(i in 1:2){
  if(i==1) pdf("results/figs/CCLE_analysis_proto.pdf")
  
  # 2 Box plots
  ################
  
  # Vectors
  HLA_vector<- as.numeric(HLA_matrix)
  mut_vector<- as.logical(mut_matrix!=0)
  mut_vector<- mut_vector[!is.na(HLA_vector)]
  HLA_vector<- HLA_vector[!is.na(HLA_vector)]
  
  # Mutations lower PHBR score? 
  cat("Obs: ",median(HLA_vector[mut_vector==T],na.rm = T),"\n") # 2.192565  
  cat("Unobs: ",median(HLA_vector[mut_vector==F],na.rm = T),"\n") # 1.544127 
  p_tmp<-  wilcox.test(HLA_vector~mut_vector)$p.value # 2.023448e-15
  
  if(i==2) svglite::svglite("results/figs/CCLE_boxplot_proto.svg")
  par(mfrow=c(2,4))
  bp<- boxplot(HLA_vector~mut_vector,outline=F,names=c("-","+"),xlab="Mutated",ylab="PHBR",frame.plot=F)
  text(1.5,max(bp$stats),labels = paste0("P=",format(p_tmp,scientific = T,digits=3)),xpd=NA)
  if(i==2) dev.off()
  
  # 3 Logreg
  ############
  source("scripts/functions/do_logreg.R")
  logreg_wp<- do_logreg(mut_matrix,HLA_matrix,"wp") # 2.24 p=9.06e-13
  logreg_wm<- do_logreg(mut_matrix,HLA_matrix,"wm") # 1.00, p= 0.92
  
  if(i==2) svglite::svglite("results/figs/CCLE_LR_proto.svg")
  par(mfrow=c(2,4))
  bp<- barplot(c(logreg_wm["OR"],logreg_wp["OR"]),ylim=c(0.95,1.30),xpd=F,ylab="Odds ratio",names.arg = c("Within mutations","Within patients"))
  arrows(bp[,1],c(logreg_wm["2.5 %"],logreg_wp["2.5 %"]),bp[,1],c(logreg_wm["97.5 %"],logreg_wp["97.5 %"]),angle = 90, code = 3,length = 0.1,xpd=NA)
  abline(h=1, lty=2)
  text(bp, rep(1.27,2),labels = signif(c(logreg_wm["p"],logreg_wp["p"]),3),xpd=NA)
  if(i==2) dev.off()
  
  # 4 PHBR ~ Frequency
  ######################
  mut_freq<- 100*colMeans(mut_matrix)
  
  # Use median PHBR per mutation
  PHBR_med<- apply(HLA_matrix,2,"median",na.rm=T)
  PHBR_freq_t<- as.data.frame(na.omit(cbind(mut_freq,PHBR=PHBR_med[names(mut_freq)])))
  
  # Plot
  if(i==2) svglite::svglite("results/figs/CCLE_mutFreq_proto.svg")
  par(mfrow=c(1,1))
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
  
  # 5 Focus on top 10
  ####################
  if(i==2) svglite::svglite("results/figs/CCLE_boxplot_top10_proto.svg")
  par(mfrow=c(1,2))
  boxplot(list(PHBR_med[1:10],PHBR_med[1:30],PHBR_med[1:100],PHBR_med),outline=F,frame.plot=F,names=c("Top 10", "Top 30", "Top 100","All"),las=2,ylab="PHBR")
  dev.off()
  
  # Stats?
  summary(PHBR_med[1:10])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.8759  1.6549  2.1889  3.4433  4.9389  9.5661 
  summary(PHBR_med)
 #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
 # 0.01931  0.56140  1.54413  2.85121  3.62629 25.11809        3 
  
} 

save(HLA_matrix,mut_matrix,file="results/data/CCLE_vp_results.RData")


