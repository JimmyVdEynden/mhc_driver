# Load data
###########
mut_matrix<- readRDS("data/mut_matrix.rds")
rownames(mut_matrix)<- substr(rownames(mut_matrix),1,12)
HLA_matrix<- readRDS("data/HLA_matrix.rds")

# 1. Plot mutation matrix
###########################

# Load cancer ids
TCGA_cancer_id <- readRDS("downloads/TCGA/clin/TCGA_cancer_id.rds")
TCGA_cancer_id<- sort(TCGA_cancer_id)

# Get mutfreq per cancer
mut_matrix_can_freq<- apply(mut_matrix[,1:20],2,function(x) tapply(x,TCGA_cancer_id[rownames(mut_matrix)],"mean",na.rm=T) )

# Only show cancer ids with data for minimal 100 samples
samples_common<- intersect(rownames(mut_matrix),rownames(HLA_matrix[!is.na(HLA_matrix[,1]),]))
samples_common<- intersect(names(TCGA_cancer_id),samples_common)
cancer_t<- table(TCGA_cancer_id[samples_common])
cancers<- names(cancer_t[cancer_t>=100])
length(cancers) # 25 cancer types
mut_matrix_can_freq<- mut_matrix_can_freq[cancers,]

# Plot
for(i in 1:2){
  if(i==1) pdf("results/figs/mut_matrix_20mut.pdf")
  if(i==2) svglite::svglite("results/figs/mut_matrix_20mut.svg")
  library(RColorBrewer)
  col_pal<- grey.colors(100,start = 0, end = 1,rev = T,gamma = 1)
  col_idx<- round(1+99*mut_matrix_can_freq/max(mut_matrix_can_freq))
  mut_matrix_cols<- apply(col_idx,2,function(x) col_pal[x])
  color2D.matplot(mut_matrix_can_freq,cellcolors = mut_matrix_cols,border = NA,axes = F,xlab = NA, ylab = NA)
  axis(1,at=(1:20)-0.5,labels = colnames(mut_matrix_can_freq)[1:20],las=2,cex.axis=0.7)
  axis(2,at=(nrow(mut_matrix_can_freq):1)-0.5,labels = rownames(mut_matrix_can_freq),las=2,cex.axis=0.7)
  dev.off()
}

# Color key
svglite::svglite("results/figs/mut_matrix_20mut_key.svg")
plot(1:100,rep(1,100),col=col_pal,pch="|",cex=5,axes=F,xlab="",ylab="")
axis(1,at = c(1,100),labels = c(0,round(100*max(mut_matrix_can_freq))))
dev.off()


# 2. Plot mut. load by sample
##################################
for(i in 1:2){
  if(i==1) pdf("results/figs/mut_load.pdf")
  if(i==2) svglite::svglite("results/figs/mut_load.svg")
  barplot(rev(rowSums(mut_matrix)),horiz = T,ylab=NA,xlab="# driver mutations",names.arg = NA)
  dev.off()
}

# 3. Plot driver mut. frequency 
##################################
for(i in 1:2){
  if(i==1) pdf("results/figs/mut_freq.pdf")
  if(i==2) svglite::svglite("results/figs/mut_freq.svg")
  barplot(100*colMeans(mut_matrix),horiz = F,xlab=NA,ylab="% mutated",names.arg = NA,border=c(rep("black",20),rep("grey",ncol(mut_matrix)-20)))
  dev.off()
}

# summary(colMeans(mut_matrix))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0004857 0.0004857 0.0006799 0.0013864 0.0010685 0.0571151 
# colMeans(mut_matrix)[colMeans(mut_matrix)>0.01]
# BRAF_V600E    IDH1_R132H  PIK3CA_E545K PIK3CA_H1047R     KRAS_G12D 
# 0.05711510    0.03846527    0.02690627    0.02515784    0.02000971 
# KRAS_G12V  PIK3CA_E542K    TP53_R175H    TP53_R273C    TP53_R248Q 
# 0.01767848    0.01699854    0.01554153    0.01282176    0.01262749 
# NRAS_Q61R     KRAS_G12C    TP53_R273H 
# 0.01155901    0.01087907    0.01058766 
# sum(colMeans(mut_matrix)>0.01) # 13
# sum(colMeans(mut_matrix)>0.001) # 186

# sort(tapply(mut_matrix[,1],TCGA_cancer_id[rownames(mut_matrix)],"mean"),decreasing=T)[1:5]
# THCA       SKCM       COAD       CHOL       LUAD 
# 0.61800000 0.43589744 0.12807882 0.05555556 0.01934236 
# 
# sort(tapply(mut_matrix[,2],TCGA_cancer_id[rownames(mut_matrix)],"mean"),decreasing=T)[1:5]
# LGG         GBM        LAML        BLCA        PRAD 
# 0.701149425 0.057500000 0.007194245 0.002433090 0.002008032 
# 
# sort(tapply(mut_matrix[,3],TCGA_cancer_id[rownames(mut_matrix)],"mean"),decreasing=T)[1:5]
# CESC       COAD       BLCA        UCS       BRCA 
# 0.13058419 0.08374384 0.07299270 0.07017544 0.06774520 

for(i in 1:2){
  if(i==1) pdf("results/figs/mut_freq_20mut.pdf")
  if(i==2) svglite::svglite("results/figs/mut_freq_20mut.svg")
  barplot(100*colMeans(mut_matrix)[1:20],horiz = F,xlab=NA,ylab="% mutated",names.arg = NA)
  abline(h=100*median(colMeans(mut_matrix)),col="red")
  dev.off()
}

