############################
# LOO analysis with logreg
#############################

# Load
mut_matrix<- readRDS("data/vp_mut_matrix.rds")
kg<- readRDS("results/data/1kg_LR.rds")

# Functions
source("scripts/functions/do_logreg.R")
library(doParallel)
do_LOO<- function(i){
  filename_tmp<- paste0("temp/LOO_1kg/",pop,"/",i,".rds")
  if(!file.exists(filename_tmp)){
    logreg_tmp<- do_logreg(mut_matrix[,-i],HLA_matrix[,-i],"wp")
    res<- c(mut=colnames(mut_matrix)[i],logreg_tmp)
    saveRDS(res,filename_tmp)
    return(res)
  }
}

for(pop in names(kg)){
  
  cat(pop," ...","\n")
  # Get HLA matrix
  HLA_matrix<- kg[[pop]][["vp_HLA1_mut_matrix"]]
  HLA_matrix<- HLA_matrix[,colnames(mut_matrix)]
  
  # Analysis
  ##########
  
  mclapply(1:ncol(mut_matrix),"do_LOO",mc.cores = 25)
  
  # Merge results
  LOO_t<- as.data.frame(matrix(NA,ncol(mut_matrix),5,dimnames = list(colnames(mut_matrix),c("mut","OR","p","ci_low","ci_high"))))
  for(i in 1:ncol(mut_matrix)){
    file_tmp<- paste0("temp/LOO_1kg/",pop,"/",i,".rds")
    if(file.exists(file_tmp)) res_tmp<- readRDS(file_tmp)
    else next
    if(is.null(res_tmp)) next
    LOO_t[i,]<- res_tmp
  }
  
  # Save
  saveRDS(LOO_t, paste0("results/data/vp_logreg_loo_1kg_",pop,".rds"))
  
}
