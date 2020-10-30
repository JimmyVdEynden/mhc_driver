# Correlate delta PHBR to mutation
################################################

# Load
load("data/PHBR_proto.RData")
mut_freq<- colSums(readRDS("data/vp_mut_matrix.rds"))
muts<- readRDS("data/driver_muts13.rds")

# dPHBR 
dPHBR_df<- as.data.frame(cbind(mut_freq,PHBR_mut=PHBR_proto_mhc1_mut,PHBR_wt=PHBR_proto_mhc1_wt))
dPHBR_df$dPHBR<- dPHBR_df$PHBR_mut-dPHBR_df$PHBR_wt

# Barplot with mut & wt for 13 muts
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_dPHBR.pdf")
  else svglite::svglite("results/figs/manuscript_dPHBR.svg",width = 15)
  # par(mfrow=c(2,1))
  barplot(t(dPHBR_df[muts,c("PHBR_wt","PHBR_mut")]),beside = T,ylab ="PHBR",las=2)
  legend("topright",c("wt peptides","mut peptides"),fill=grey.colors(2), bty="n")
  dev.off()
}

# How does BRAF relate to other subst? ECDF
for(i in 1:2){
  if(i==1) pdf("results/figs/manuscript_dPHBR_ecdf.pdf")
  else svglite::svglite("results/figs/manuscript_dPHBR_ecdf.svg")
  par(mfrow=c(1,2))
  dPHBR_ecdf<- ecdf(dPHBR_df$dPHBR)
  plot(dPHBR_ecdf,frame.plot=F,xlab="dPHBR",ylab = "Cumulative frequency")
  # dPHBR_ecdf(sort(dPHBR_df$dPHBR))
  for(mut in muts){
    mut_co<- c(dPHBR_df[mut,"dPHBR"],dPHBR_ecdf(dPHBR_df[mut,"dPHBR"]))
    points(mut_co[1],mut_co[2],pch=16,col="red")
    text(mut_co[1],mut_co[2],gsub("_"," ",mut),cex=0.5,ad=c(0,1),col="red")
  }
  dev.off()
}
cat(muts[1],c(dPHBR_df[muts[1],"dPHBR"],dPHBR_ecdf(dPHBR_df[muts[1],"dPHBR"])))
cat(muts[11],c(dPHBR_df[muts[11],"dPHBR"],dPHBR_ecdf(dPHBR_df[muts[11],"dPHBR"])))
# BRAF V600E:  7.655281 0.9956012
# TP53_Y220C 3.510428 0.9442815
