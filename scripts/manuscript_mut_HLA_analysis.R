# Load
HLA_matrix<- readRDS("data/HLA_matrix.rds")
mut_matrix<- readRDS("data/mut_matrix.rds")
TCGA_cancer_id<- readRDS("downloads/TCGA/clin/TCGA_cancer_id.rds")

# Get score vector
HLA_vector<- as.numeric(HLA_matrix)
mut_vector<- as.logical(mut_matrix!=0)
mut_vector<- mut_vector[!is.na(HLA_vector)]
HLA_vector<- HLA_vector[!is.na(HLA_vector)]

# Mutations lower PHBR score? 
median(HLA_vector[mut_vector==T],na.rm = T) # 2.08
median(HLA_vector[mut_vector==F],na.rm = T) # 1.58
p_tmp<-  wilcox.test(HLA_vector~mut_vector)$p.value # e-95

# Randomize patients: same effect!!!
mut_matrix_rand_pt<- mut_matrix
for(i in 1:ncol(mut_matrix)){
  mut_matrix_rand_pt[,i]<- sample(mut_matrix_rand_pt[,i])
}

# Randomize mutations: effect gone!!!
mut_matrix_rand_mut<- mut_matrix
for(i in 1:nrow(mut_matrix)){
  mut_matrix_rand_mut[i,]<- sample(mut_matrix_rand_mut[i,])
}

# Bring together in 1 plot
##########################

mut_matrix_ls<- list(mut_matrix,mut_matrix_rand_mut,mut_matrix_rand_pt)
# pdf("results/figs/mut_HLA_all_boxplot.pdf")
svglite::svglite("results/figs/mut_HLA_all_boxplot.svg")
par(mfrow=c(2,4))
for(i in 1:3){
  cat(i," ")
  mut_matrix<- mut_matrix_ls[[i]]
  HLA_vector<- as.numeric(HLA_matrix)
  mut_vector<- as.logical(mut_matrix!=0)
  mut_vector<- mut_vector[!is.na(HLA_vector)]
  HLA_vector<- HLA_vector[!is.na(HLA_vector)]
  p_tmp<-  wilcox.test(HLA_vector~mut_vector)$p.value # e-95
  if(i==1) title<- "No randomization"
  if(i==2) title<- "Mutations randomized"
  if(i==3) title<- "Patients randomized"
  boxplot(HLA_vector~mut_vector,outline=F,names=c("-","+"),xlab="Mutated?",ylab="PHBR",ylim=c(0,10),frame.plot=F,main=title)
  text(1.5,10,labels = paste0("P=",format(p_tmp,scientific = T,digits=3)),xpd=NA)
  cat(title,"\n")
  cat("Obs. ",median(HLA_vector[mut_vector==T],na.rm = T),"\n") # 2.08
  cat("Unobs. ",median(HLA_vector[mut_vector==F],na.rm = T),"\n") # 1.58
  cat("P ",p_tmp,"\n","\n")
}
dev.off()
