#############################################
# manuscript_logreg_aa_HLABinder.R
#############################################

# Load
GPPM_rand <- readRDS("data/GPPM_rand.rds")
load("../immunoediting_2019/data/aa_classes.RData")

# Create matrix with n aa
aa_matrix<- matrix(0,length(GPPM_rand),20,dimnames = list(GPPM_rand$mut,aa_all))
pep21<- GPPM_rand$pep21_mut
for(i in 1:nrow(aa_matrix)){
  cat(i," ")
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

# Linear regression on log(PHBR) 
#################################

lm_t<- data.frame(row.names = aa_all)
lm_t$p<- NA 
lm_t$b0<- NA 
lm_t$b1<- NA 
for(i in 1:2){
  if(i==1) svglite::svglite(paste0("results/figs/manuscript_lm_aa_HLABinder.svg"))
  else pdf(paste0("results/figs/manuscript_lm_aa_HLABinder.pdf"))
  for(i in 1:length(aa_all)){
    aa<- aa_all[i]
    lty_aa<- 1
    if(aa%in%aa_hp) col_aa<- "green"
    else if(aa%in%aa_polar) col_aa<- "blue"
    else col_aa<- "red"
    
    if(i==1){
      plot(NA,xlim=c(0,100),ylim=c(-2,2),frame.plot=F,axes=F,xlab="% in peptides",ylab="PHBR")
      axis(2,at=seq(-2,2,1),c(0.01,0.1,1,10,100),las=2)
      axis(1)
      # axis(1,at=seq(0,366,366/5),labels = seq(0,100,20))
    }
    lm_tmp<- lm(log10(GPPM_rand$PHBR_mut)~as.numeric(aa_matrix[,aa]/3.66))
    lm_tmp_p<- coefficients(summary(lm_tmp))[2,"Pr(>|t|)"]
    lm_tmp_b0<- coefficients(lm_tmp)[1]
    lm_tmp_b1<- coefficients(lm_tmp)[2]
    lm_t[i,c("p","b0","b1")]<- c(lm_tmp_p,lm_tmp_b0,lm_tmp_b1)
    if(lm_tmp_p>0.0025) lty_aa<- 2 # Bonferroni at 0.05
    abline(lm_tmp,lwd=3,col=col_aa,lty=lty_aa)
    text(x = 100,y = lm_tmp_b0 + lm_tmp_b1*100,labels = aa,xpd=NA,col=col_aa)
  }
  dev.off()
}
lm_t$q<- p.adjust(lm_t$p,"fdr")
lm_t
#         p          b0           b1             q
# A  3.879050e-02 0.173992704 -0.001421706  4.310055e-02
# C  2.157110e-43 0.127284280  0.014881520  4.314220e-43
# D  1.049395e-99 0.078865126  0.017828429  2.998273e-99
# E 2.420412e-274 0.006756543  0.023240178 4.840824e-273
# F 2.159902e-163 0.258725266 -0.025435363 1.079951e-162
# G 6.725121e-140 0.056505880  0.016977305 2.690048e-139
# H  4.747845e-03 0.173117280 -0.002919739  5.585700e-03
# I  1.441085e-07 0.184692742 -0.004323618  2.217053e-07
# K  2.153556e-75 0.089753080  0.013291687  4.785680e-75
# L 2.274503e-246 0.355489051 -0.020055206 2.274503e-245
# M 9.258557e-133 0.229407633 -0.029557481 3.086186e-132
# N  2.465653e-21 0.131144984  0.008820864  4.109422e-21
# P  2.534917e-05 0.181398724 -0.002829719  3.168646e-05
# Q  2.294835e-39 0.114888765  0.010685153  4.172428e-39
# R  5.551966e-98 0.254067360 -0.015173000  1.387992e-97
# S  5.790484e-06 0.141997126  0.002703112  7.720645e-06
# T  1.025056e-01 0.157510508  0.001299669  1.079006e-01
# V  4.036148e-07 0.187626867 -0.003917956  5.765926e-07
# W  1.313222e-01 0.167781233 -0.002264035  1.313222e-01
# Y 3.320013e-229 0.262146889 -0.033819313 2.213342e-228
