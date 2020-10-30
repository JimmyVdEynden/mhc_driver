# Load
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")
load("data/PHBR_proto.RData")

# 1) Mut mhc1 peptides
vp_HLA1_mut_matrix<- vp_mut_matrix
vp_HLA1_mut_matrix[,]<- NA
common_muts<- intersect(names(PHBR_proto_mhc1_mut),colnames(vp_HLA1_mut_matrix))
for(i in 1:nrow(vp_HLA1_mut_matrix)){
  vp_HLA1_mut_matrix[i,common_muts]<- PHBR_proto_mhc1_mut[common_muts]
}

# 2) WT mhc1 peptides
vp_HLA1_wt_matrix<- vp_mut_matrix
vp_HLA1_wt_matrix[,]<- NA
common_muts<- intersect(names(PHBR_proto_mhc1_wt),colnames(vp_HLA1_wt_matrix))
for(i in 1:nrow(vp_HLA1_wt_matrix)){
  vp_HLA1_wt_matrix[i,common_muts]<- PHBR_proto_mhc1_wt[common_muts]
}

# 3) Mut mhc2 peptides
vp_HLA2_mut_matrix<- vp_mut_matrix
vp_HLA2_mut_matrix[,]<- NA
common_muts<- intersect(names(PHBR_proto_mhc2_mut),colnames(vp_HLA2_mut_matrix))
for(i in 1:nrow(vp_HLA2_mut_matrix)){
  vp_HLA2_mut_matrix[i,common_muts]<- PHBR_proto_mhc2_mut[common_muts]
}

# Save
save(vp_HLA1_mut_matrix,vp_HLA1_wt_matrix,vp_HLA2_mut_matrix,file="data/vp_HLA_matrix.RData")

