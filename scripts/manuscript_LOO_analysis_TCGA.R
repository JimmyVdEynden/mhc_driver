############################
# LOO analysis with logreg
#############################

# Load
mut_matrix<- readRDS("data/mut_matrix.rds")
HLA_matrix<- readRDS("data/HLA_matrix.rds")

# Functions
source("scripts/functions/do_logreg.R")

# Analysis
##########

library(doParallel)

# Alternative
do_LOO<- function(i){
  filename_tmp<- paste0("temp/LOO_TCGA//",i,".rds")
  if(!file.exists(filename_tmp)){
    logreg_tmp<- do_logreg(mut_matrix[,-i],HLA_matrix[,-i],"wp")
    res<- c(mut=colnames(mut_matrix)[i],logreg_tmp)
    saveRDS(res,filename_tmp)
    return(res)
  }
}
mclapply(1:ncol(mut_matrix),"do_LOO",mc.cores = 25)

# Merge results
LOO_t<- as.data.frame(matrix(NA,ncol(mut_matrix),5,dimnames = list(colnames(mut_matrix),c("mut","OR","p","ci_low","ci_high"))))
for(i in 1:ncol(mut_matrix)){
  cat(i," ")
  res_tmp<- readRDS(paste0("temp/LOO_TCGA//",i,".rds"))
  if(is.null(res_tmp)) next
  LOO_t[i,]<- res_tmp
}

# Save
saveRDS(LOO_t, "results/data/vp_logreg_loo_TCGA.rds")

