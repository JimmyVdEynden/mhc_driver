# Select genes based on >1% effect
##################################

# Baseline OR
baseline<- as.data.frame(readRDS("results/data/vp_logreg_wp.rds"))["mut",]

# LOO
loo<- readRDS("results/data/vp_logreg_loo.rds")
loo<- loo[order(loo$OR),]

# Calculate % contribution based on LOO
loo$rel_effect<- 1-(loo$OR - 1)/(baseline$OR-1)

# Get & save mutations
muts<- loo$mutation[loo$rel_effect>0.01]
saveRDS(muts, "data/driver_muts13.rds")
