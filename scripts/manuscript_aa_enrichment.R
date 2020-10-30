# Load aa classes
load("../immunoediting_2019/data/aa_classes.RData")

# Load driver muts
mut_matrix<- readRDS("data/vp_mut_matrix.rds")
muts<- colnames(mut_matrix)
mut_freq<- colMeans(mut_matrix)

# Get peptides
pep21<- rep(NA, length(muts))
names(pep21)<- muts
for(i in 1:length(pep21)){
  mut_tmp<- names(pep21)[i]
  file_tmp<- paste0("temp/netMHCPan40/input/",mut_tmp,".fasta")
  if(!file.exists(file_tmp)) next
  pep21[i]<- as.character(read.table(file = file_tmp,stringsAsFactors = F,header=T))
}

# Create aa matrix
pep21[nchar(pep21)!=21]<- NA
aa_matrix<- matrix(0,length(pep21),20,dimnames = list(names(pep21),aa_all))
substr(pep21,11,11)<- "X" # leave out mutated aa to focus on context only
names(pep21)<- rownames(aa_matrix)
for(i in 1:nrow(aa_matrix)){
  pep_tmp<- NULL
  # 8-mers
  for(j in 4:11){
    pep_tmp<- c(pep_tmp,substr(pep21[i],j,j+7))
  }
  # 9-mers
  for(j in 3:11){
    pep_tmp<- c(pep_tmp,substr(pep21[i],j,j+8))
  }
  # 10-mers
  for(j in 2:11){
    pep_tmp<- c(pep_tmp,substr(pep21[i],j,j+9))
  }
  # 11-mers
  for(j in 1:11){
    pep_tmp<- c(pep_tmp,substr(pep21[i],j,j+10))
  }
  pep_tmp<- table(unlist(strsplit(pep_tmp,"")))
  pep_tmp<- pep_tmp[names(pep_tmp)!="X"]
  aa_matrix[i,names(pep_tmp)]<- pep_tmp 
}

# General correlation for aa?
aa_cor_t<- data.frame(row.names = aa_all)
aa_cor_t$p<- NA
aa_cor_t$r<- NA
for(aa in aa_all){
  freq_aa_t<- data.frame(freq=log10(mut_freq[rownames(aa_matrix)]),aa=aa_matrix[,aa])
  freq_aa_t<- freq_aa_t[!is.infinite(freq_aa_t$freq),]
  cor_tmp<- cor.test(freq_aa_t$freq,freq_aa_t$aa,method = "spearman")
  aa_cor_t[aa,c("p","r")]<- c(cor_tmp$p.value,cor_tmp$estimate)
}
aa_cor_t$q<- p.adjust(aa_cor_t$p,"fdr") 
# Only sign for Gly: r=0.16, P=1.461746e-05, q=0.000293

# Plot
for(i in 1:2){
  if(i==1) svglite::svglite(paste0("results/figs/manuscript_aa_enrichment_G_correlation.svg"))
  else pdf(paste0("results/figs/manuscript_aa_enrichment_G_correlation.pdf"))
  aa<- "G"
  freq_aa_t<- data.frame(freq=log10(mut_freq[rownames(aa_matrix)]),aa=aa_matrix[,aa]/3.66) # /366 to normalize to max number
  freq_aa_t<- freq_aa_t[!is.infinite(freq_aa_t$freq),]
  plot(freq_aa_t$freq,freq_aa_t$aa,pch=16,ylim=c(0,25),frame.plot=F,axes=F,xlab="Mutation frequency (%)",ylab="Gly residues in peptide (%)")
  axis(2)
  axis(1,at=c(-4,-3,-2,-1),labels = c(0.01,0.1,1,10))
  abline(lm(freq_aa_t$aa~freq_aa_t$freq),lwd=3,col="red")
  dev.off()
}

# Enrichment for each mutation?
p_matrix<- aa_matrix
p_matrix[,]<- NA
q_matrix<- p_matrix
OR_matrix<- p_matrix
for(i in 1:nrow(aa_matrix)){
  cat(i," ")
  for(aa in colnames(aa_matrix)){
    idxAA<- which(colnames(aa_matrix)==aa)
    aa_t<- cbind(c(sum(aa_matrix[i,idxAA]),sum(aa_matrix[-i,idxAA])),c(sum(aa_matrix[i,-idxAA]),sum(aa_matrix[-i,-idxAA])))
    fish_tmp<- fisher.test(aa_t,alternative = "greater")
    p_matrix[i,idxAA]<- fish_tmp$p.value
    OR_matrix[i,idxAA]<- fish_tmp$estimate
  }
  q_matrix[i,]<- p.adjust(p_matrix[i,],"fdr")
}
# View(t(p_matrix))
# View(t(OR_matrix))

OR_matrix[c("KRAS_G12D","KRAS_G12V","KRAS_G12C","KRAS_G13D"),c("G")]
p_matrix[c("KRAS_G12D","KRAS_G12V","KRAS_G12C","KRAS_G13D"),c("G")]

OR_matrix[c("TP53_H179R","TP53_R175H"),c("C","H")]
p_matrix[c("TP53_H179R","TP53_R175H"),c("C","H")]

# Heatmap?
muts13<- readRDS("data/driver_muts13.rds")
# data_hm<- log(OR_matrix[muts13,]+0.01)
data_hm<- OR_matrix[muts13,]
data_hm[q_matrix[muts13,]>0.05]<- 0
for(i in 1:2){
  if(i==1) svglite::svglite(paste0("results/figs/manuscript_aa_enrichment_hm.svg"))
  else pdf("results/figs/manuscript_aa_enrichment_hm.pdf")
  hm<- heatmap.2(data_hm,col=bluered(75)[38:75],symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none',distfun=function(x) dist(x,method = "manhattan"))
  dev.off()
}

# Key
svglite::svglite("results/figs/manuscript_aa_enrichment_hm_key.svg")
par(mfrow=c(3,3))
plot(1:38,rep(1,38),pch="|",col=bluered(75)[38:75],cex=5,axes=F,xlab=NA,ylab=NA)
axis(1,at = c(1,4*38/max(data_hm,na.rm=T),8*38/max(data_hm,na.rm=T)),labels = c(0,4,8))
dev.off()

cat(paste(muts13[rev(hm$rowInd)],pep21[muts13][rev(hm$rowInd)]),sep="\n")
# KRAS_G12C TEYKLVVVGAXGVGKSALTIQ
# KRAS_G12D TEYKLVVVGAXGVGKSALTIQ
# KRAS_G12V TEYKLVVVGAXGVGKSALTIQ
# KRAS_G13D EYKLVVVGAGXVGKSALTIQL
# TP53_Y220C NTFRHSVVVPXEPPEVGSDCT
# PIK3CA_E545K STRDPLSEITXQEKDFLWSHR
# PIK3CA_E542K KAISTRDPLSXITEQEKDFLW
# TP53_R282W VRVCACPGRDXRTEEENLRKK
# BRAF_V600E VKIGDFGLATXKSRWSGSHQF
# IDH1_R132H SGWVKPIIIGXHAYGDQYRAT
# TP53_R175H QSQHMTEVVRXCPHHERCSDS
# TP53_H179R MTEVVRRCPHXERCSDSDGLA
# FBXW7_R465C HTLYGHTSTVXCMHLHEKRVV