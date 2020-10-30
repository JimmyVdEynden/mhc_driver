#######################################################
# Get within-atient odds ratios for wt and mut VP data
#######################################################

# logreg function
source("scripts/functions/do_logreg.R")

# Load
vp_mut_matrix<- readRDS("data/vp_mut_matrix.rds")
load("data/vp_HLA_matrix.RData")

# MUT
######

cat("MUT ...", "\n")

# Select matrices
vp_HLA_matrix_sel<- vp_HLA1_mut_matrix
vp_mut_matrix_sel<- vp_mut_matrix[,colnames(vp_HLA_matrix_sel)]

# Logreg
logreg_mut<- do_logreg(vp_mut_matrix_sel,vp_HLA_matrix_sel,"wp")

# WT
####

cat("WT ...", "\n")

# Select matrices
vp_HLA_matrix_sel<- vp_HLA1_wt_matrix
vp_mut_matrix_sel<- vp_mut_matrix[,colnames(vp_HLA_matrix_sel)]

# Logreg
logreg_wt<- do_logreg(vp_mut_matrix_sel,vp_HLA_matrix_sel,"wp")

# Save
######
vp_logreg_wp<- rbind(mut=logreg_mut,wt=logreg_wt)
saveRDS(vp_logreg_wp,"results/data/vp_logreg_wp.rds")