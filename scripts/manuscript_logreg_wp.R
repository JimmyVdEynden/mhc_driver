# # Within patient analysis per cancer type
############################################

# Functions
source("scripts/functions/do_logreg.R")

# Data
data_per_cancer_type<- readRDS("data/data_per_cancer_type.rds")

# Restrict to cancer types with minimal 100 samples in both datasets
cancers<- names(data_per_cancer_type)[sapply(data_per_cancer_type,function(x) nrow(x$mutations)>=100)]

# Analysis
loo_wp_t<- as.data.frame(matrix(NA,length(cancers),4,dimnames = list(cancers,c("OR","p","ci_low","ci_high"))))

for(cancer in cancers){
  cat(cancer,"\n")
  mut_matrix_sel<- data_per_cancer_type[[cancer]]$mutations
  HLA_matrix_sel<- data_per_cancer_type[[cancer]]$affinities
  loo_wp_t[cancer,]<- do_logreg(mut_matrix_sel,HLA_matrix_sel,"wp")
}

# Save
saveRDS(loo_wp_t, "results/data/logreg_wp.rds")
