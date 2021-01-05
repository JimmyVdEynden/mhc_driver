#######################################################
# Get within-atient odds ratios for wt and mut VP data
#######################################################

# logreg function
source("scripts/functions/do_logreg.R")

# Load
mut_matrix<- readRDS("data/mut_matrix.rds")
HLA_matrix<- readRDS("data/HLA_matrix.rds")
HLA_wt_matrix<- readRDS("data/HLA_wt_matrix.rds")

# MUT
######

cat("MUT ...", "\n")

# Logreg
logreg_mut<- do_logreg(mut_matrix,HLA_matrix,"wp")

# WT
####

cat("WT ...", "\n")

# Logreg
logreg_wt<- do_logreg(mut_matrix,HLA_wt_matrix,"wp")

# Save
######
logreg_wp<- rbind(mut=logreg_mut,wt=logreg_wt)
saveRDS(logreg_wp,"results/data/TCGA_logreg_wp.rds")

