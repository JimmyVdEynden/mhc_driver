# Load data
driver_overview<- readRDS(file="data/driver_overview_cons.rds")
muts<- readRDS("data/driver_muts13.rds")

# Conditions
driver_overview$cond<- NA
driver_overview$cond[driver_overview$isObserved==F]<- "Unobserved"
# driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==T]<- "Rec_D"
driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==T&rownames(driver_overview)%in%muts]<- "Rec_D_HA"
driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==T&!rownames(driver_overview)%in%muts]<- "Rec_D_noHA"
driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==F]<- "noRec_D"
driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==F&driver_overview$isRecurr==T]<- "Rec_noD"
driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==F&driver_overview$isRecurr==F]<- "noRec_noD"
driver_overview$cond<- factor(driver_overview$cond,levels=c("Unobserved","noRec_noD","Rec_noD","noRec_D","Rec_D_noHA","Rec_D_HA"))

# driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==T&!rownames(driver_overview)%in%muts&driver_overview$isTTCG]<- "Rec_D_noHA_TTCG"
# driver_overview$cond[driver_overview$isObserved==T&driver_overview$isCGC==T&driver_overview$isRecurr==T&!rownames(driver_overview)%in%muts&!driver_overview$isTTCG]<- "Rec_D_noHA_noTTCG"

# Plot
svglite::svglite("results/figs/manuscript_cons.svg",height = 10)
par(mfcol=c(3,2))
# PHBR
boxplot(driver_overview$PHBR~driver_overview$cond,outline=F,xlab="",ylab="PHBR",frame.plot=F)
wilcox.test(driver_overview$PHBR[driver_overview$cond=="Rec_D_HA"],driver_overview$PHBR[driver_overview$cond=="Rec_D_noHA"]) # P= 3.279e-05
# Compare phastcons
boxplot(driver_overview$phastCons100way~driver_overview$cond,outline=F,xlab="",ylab="phastCons100way",frame.plot=F)
wilcox.test(driver_overview$phastCons100way[driver_overview$cond=="Rec_D_HA"],driver_overview$phastCons100way[driver_overview$cond=="Rec_D_noHA"]) # P= 0.17
# Compare phyloP
boxplot(driver_overview$phyloP100way~driver_overview$cond,outline=F,xlab="",ylab="phyloP100way",frame.plot=F)
wilcox.test(driver_overview$phyloP100way[driver_overview$cond=="Rec_D_HA"],driver_overview$phyloP100way[driver_overview$cond=="Rec_D_noHA"]) # P= 0.21; 0.02318 if +/-10 regions
dev.off()

# Numbers
tapply(driver_overview$phastCons100way,driver_overview$cond,"summary")
# $Rec_D_noHA
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.8095  0.9143  0.8453  0.9714  1.0000 
# 
# $Rec_D_HA
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5571  0.7571  0.9524  0.8930  1.0000  1.0000 

tapply(driver_overview$phyloP100way,driver_overview$cond,"summary")
# $Rec_D_noHA
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.452   3.762   4.881   4.515   5.643   7.476 
# 
# $Rec_D_HA
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.381   4.381   5.976   5.454   6.429   6.929