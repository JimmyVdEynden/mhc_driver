############################
# LOO analysis with logreg
#############################

# library(tidyverse)

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

LOO_t<- as.data.frame(matrix(NA,ncol(vp_mut_matrix_sel),4,dimnames = list(colnames(vp_mut_matrix_sel),c("OR","p","ci_low","ci_high"))))

# for(i in 1:nrow(LOO_t)){
#   loo_wp_t[i,]<- do_logreg_wp(vp_mut_matrix_sel[,-i],vp_HLA_matrix_sel[,-i])
# }

library(foreach)
library(doParallel)
cl <- makeCluster(50) 
registerDoParallel(cl)
# foreach(i=1:nrow(LOO_t)) %dopar% {
  foreach(i=1:10) %dopar% {
    cat(i," ")
    LOO_t[i,]<- do_logreg(vp_mut_matrix_sel[,-i],vp_HLA_matrix_sel[,-i],"wp")
  i
}
stopCluster(cl)

# Save
saveRDS(loo_t, "results/data/vp_logreg_loo.rds")

