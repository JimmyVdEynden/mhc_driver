# OR after leaving 13 selected mutations out?
#################################################

# logreg function
source("scripts/functions/do_logreg_wp.R")

# Load muts
muts<- readRDS("data/driver_muts13.rds")

# Load
vp_mut_matrix<- readRDS("data/vp_mut_matrix.rds")
load("data/vp_HLA_matrix.RData")

# Remove muts
idx_excl<- which(colnames(vp_HLA1_mut_matrix)%in%muts)
vp_HLA_matrix<- vp_HLA1_mut_matrix[,-idx_excl]
vp_mut_matrix<- vp_mut_matrix[,colnames(vp_HLA_matrix)]

# Logreg
logreg_excl_mut13<- do_logreg_wp(vp_mut_matrix,vp_HLA_matrix)

# Save
saveRDS(logreg_excl_mut13,"results/data/vp_logreg_loo_mut13.rds")
