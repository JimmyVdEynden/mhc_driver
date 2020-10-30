do_logreg<- function(mut_matrix,HLA_matrix,model=c("wp","wm")){
  
  # Remove NA columns
  idx_NA_col<- which(is.na(HLA_matrix[1,]))
  if(length(idx_NA_col)!=0){
    mut_matrix_tmp<- mut_matrix[,-idx_NA_col]
    HLA_matrix_tmp<- HLA_matrix[,-idx_NA_col]
  }
  else{
    mut_matrix_tmp<- mut_matrix
    HLA_matrix_tmp<- HLA_matrix
  }
  
  # get vectors
  mut_vector<- as.logical(mut_matrix_tmp != 0)
  HLA_vector<- as.numeric(HLA_matrix_tmp)
  
  # Calculate
  if(model=="wp") MR_effect<- rep(rowMeans(mut_matrix_tmp),ncol(mut_matrix_tmp))
  if(model=="wm") MR_effect<- rep(colMeans(mut_matrix_tmp),each=nrow(mut_matrix_tmp))
  glmer_pt <- lme4::glmer(mut_vector~log(HLA_vector) + (1|MR_effect),family=binomial)
  OR_tmp<- exp(summary(glmer_pt)$coefficients["log(HLA_vector)","Estimate"])
  p_tmp<- summary(glmer_pt)$coefficients["log(HLA_vector)","Pr(>|z|)"]
  
  ci<- confint(glmer_pt, level = 0.95, method = "Wald")
  ci<- exp(ci[grep("log",rownames(ci)),])
  res<- c(OR=OR_tmp, p=p_tmp, ci)
  return(res)
  
}
