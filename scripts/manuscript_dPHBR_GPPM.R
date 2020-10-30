################################################################
# Correlate delta PHBR to substitution for random substitutions
################################################################

# Load
GPPM_rand<- readRDS(file="data/GPPM_rand.rds")

# Create substitution matrix with median PHBR
load("../immunoediting_2019/data/aa_classes.RData")
aa_dPHBR_matrix<- matrix(NA,length(aa_all),length(aa_all),dimnames = list(aa_all,aa_all))
for(i in 1:nrow(aa_dPHBR_matrix)){
  cat(i," ")
  ref_aa_tmp<- rownames(aa_dPHBR_matrix)[i]
  GPPM_temp<- GPPM_rand[GPPM_rand$ref_aa==ref_aa_tmp]
  alt_aa_PHBRm_tmp<- tapply(GPPM_temp$mhc1_mut_mean_aff_rank,factor(GPPM_temp$alt_aa,levels=aa_all),"median",na.rm=T)
  alt_aa_PHBRwt_tmp<- tapply(GPPM_temp$mhc1_wt_mean_aff_rank,factor(GPPM_temp$alt_aa,levels=aa_all),"median",na.rm=T)
  aa_dPHBR_matrix[ref_aa_tmp,]<- alt_aa_PHBRm_tmp-alt_aa_PHBRwt_tmp
}

# Heatmap
##########

# Row order on medians
ref_med<- rowMedians(aa_dPHBR_matrix,na.rm=T)
names(ref_med)<- colnames(aa_dPHBR_matrix)
# sort(ref_med)

# Col order on medians
alt_med<- rowMedians(t(aa_dPHBR_matrix),na.rm=T)
names(alt_med)<- rownames(aa_dPHBR_matrix)
# sort(alt_med)

# Plot hm
svglite::svglite("results/figs/manuscript_dPHBR_GPPM_heatmap.svg")
data_hm<- aa_dPHBR_matrix
data_hm[is.na(data_hm)]<- 0
data_hm<- data_hm[order(ref_med),order(alt_med)]
heatmap.2(data_hm,col=bluered(75),symkey=TRUE, key=F, keysize=1, trace="none",scale="none",density.info='none',Rowv=NULL,Colv=NULL)
dev.off()

# Key
svglite::svglite("results/figs/manuscript_dPHBR_GPPM_heatmap_key.svg")
par(mfrow=c(3,3))
plot(1:75,rep(1,75),pch="|",col=bluered(75),cex=5,axes=F,xlab=NA,ylab=NA)
axis(1,at = c(37.5-4*37.5/max(aa_dPHBR_matrix,na.rm=T),37.5,37.5+4*37.5/max(aa_dPHBR_matrix,na.rm=T)),labels = c(-4,0,4))
dev.off()

# Boxplots 
svglite::svglite("results/figs/manuscript_dPHBR_GPPM_bp.svg")
par(mfrow=c(1,2))
# ref
data_hm<- aa_dPHBR_matrix
data_hm<- data_hm[order(ref_med,decreasing = T),order(alt_med,decreasing = T)]
boxcols<- bluered(75)
# col_idx<- round(75*(alt_med-min(alt_med))/(max(alt_med)-min(alt_med)))
col_idx<- round(75*(alt_med-min(aa_dPHBR_matrix,na.rm=T))/(max(aa_dPHBR_matrix,na.rm=T)-min(aa_dPHBR_matrix,na.rm=T)))
col_idx[col_idx==0]<- 1
boxcols<- boxcols[col_idx]
boxcols<- boxcols[order(alt_med,decreasing = T)] 
# boxplot(data_hm,horizontal=F, outline=F,staplewex=0,lty=1,frame.plot=F,xlab="dPHBR",las=2,col=boxcols)
boxplot(data_hm,horizontal=T, outline=F,staplewex=0,lty=1,frame.plot=F,xlab="dPHBR",las=2,col=boxcols,main="Alt aa")
# alt
boxcols<- bluered(75)
# col_idx<- round(75*(ref_med-min(ref_med))/(max(ref_med)-min(ref_med)))
col_idx<- round(75*(ref_med-min(aa_dPHBR_matrix,na.rm=T))/(max(aa_dPHBR_matrix,na.rm=T)-min(aa_dPHBR_matrix,na.rm=T)))
col_idx[col_idx==0]<- 1
boxcols<- boxcols[col_idx]
boxcols<- boxcols[order(ref_med,decreasing = T)] 
boxplot(t(data_hm),horizontal=T, outline=F,staplewex=0,lty=1,frame.plot=F,xlab="dPHBR",las=2,col=boxcols,main="Ref aa")
dev.off()

# PP2 
svglite::svglite("results/figs/manuscript_dPHBR_GPPM_PP2.svg")
par(mfrow=c(1,4))
# PP2 ref
vioplot::vioplot(GPPM_rand$PP2~factor(GPPM_rand$ref_aa,levels=rev(rownames(data_hm))),frame.plot=F,ylab=NA,xlab="PP2",horizontal=T,las=2,col="grey85",colMed='black',rectCol=NA,lineCol=NA)
# PP2 alt
vioplot::vioplot(GPPM_rand$PP2~factor(GPPM_rand$alt_aa,levels=rev(colnames(data_hm))),frame.plot=F,ylab=NA,xlab="PP2",horizontal=T,las=2,col="grey85",colMed='black',rectCol=NA,lineCol=NA)
dev.off()

# Some numbers
cat(aa_dPHBR_matrix["V","E"]) # 4.482721
cat(alt_med["E"]) # 1.76
cat(alt_med["D"]) # 2.41
cat(ref_med["Y"]) # 2.17
tapply(GPPM_rand$PP2,factor(GPPM_rand$ref_aa,levels=rev(rownames(data_hm))),"summary")
