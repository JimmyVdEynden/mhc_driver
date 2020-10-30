# Create virtual pt data
########################

# Mut matrix
############

# Driver mut probability (pan cancer)
mut_matrix<- readRDS("data/mut_matrix.rds")
Pmut<- colMeans(mut_matrix)

# Create virtual pts
vp_mut_matrix<- matrix(NA,nrow(mut_matrix),ncol(mut_matrix))
colnames(vp_mut_matrix)<- colnames(mut_matrix)
for(i in 1:nrow(vp_mut_matrix)){
  cat(i," ")
  for(j in 1:ncol(vp_mut_matrix)){
    vp_mut_matrix[i,j]<- sample(x = c(0,1),size = 1,replace = T, prob = c(1-Pmut[j],Pmut[j]))
  }
}

dim(vp_mut_matrix)
# [1]] 10295   688

# How many BRAF?
sum(vp_mut_matrix[,"BRAF_V600E"])

# # Check
# barplot(colSums(vp_mut_matrix)[1:50])
# 
# Save
saveRDS(vp_mut_matrix,file="data/vp_mut_matrix.rds")

