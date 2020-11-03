############################
# LOO analysis with logreg
#############################

# Load
vp_mut_matrix<- readRDS("data/vp_mut_matrix.rds")
load("data/vp_HLA_matrix.RData")

# Functions
source("scripts/functions/do_logreg.R")

# Select matrices
vp_HLA_matrix_sel<- vp_HLA1_mut_matrix
vp_mut_matrix_sel<- vp_mut_matrix[,colnames(vp_HLA_matrix_sel)]

# Analysis
##########

# library(foreach)
library(doParallel)
# cl <- makeCluster(25,outfile="log/LOO_log.txt") 
# registerDoParallel(cl)
# LOO_t<- foreach(i=1:ncol(vp_mut_matrix_sel)) %dopar% {
#   filename_tmp<- paste0("temp/LOO/",i,".rds")
#   if(file.exists(filename_tmp)) next
#   logreg_tmp<- do_logreg(vp_mut_matrix_sel[,-i],vp_HLA_matrix_sel[,-i],"wp")
#   res<- c(mut=colnames(vp_mut_matrix_sel)[i],logreg_tmp)
#   saveRDS(res,filename_tmp)
#   return(res)
# }
# stopCluster(cl)

# Alternative
do_LOO<- function(i){
  filename_tmp<- paste0("temp/LOO/",i,".rds")
  if(!file.exists(filename_tmp)){
    logreg_tmp<- do_logreg(vp_mut_matrix_sel[,-i],vp_HLA_matrix_sel[,-i],"wp")
    res<- c(mut=colnames(vp_mut_matrix_sel)[i],logreg_tmp)
    saveRDS(res,filename_tmp)
    return(res)
  }
}
mclapply(1:ncol(vp_mut_matrix_sel),"do_LOO",mc.cores = 25)

# Merge results
LOO_t<- as.data.frame(matrix(NA,ncol(vp_mut_matrix_sel),5,dimnames = list(colnames(vp_mut_matrix_sel),c("mut","OR","p","ci_low","ci_high"))))
for(i in 1:ncol(vp_mut_matrix_sel)){
  res_tmp<- readRDS(paste0("temp/LOO/",i,".rds"))
  if(is.null(res_tmp)) next
  LOO_t[i,]<- res_tmp
}

# Save
saveRDS(LOO_t, "results/data/vp_logreg_loo.rds")

