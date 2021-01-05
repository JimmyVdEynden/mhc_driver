source("scripts/functions/do_logreg.R")

# Get drive rmutation
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))# 
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")

for(m in 1:length(driver_muts)){
  cat(m," ")
  mut_tmp<- driver_muts[m]
  gene_mhc1_filename<- paste0("temp/netMHCPan40/output/",mut_tmp,1,'.xls') 
  if(!file.exists(gene_mhc1_filename)) next
  for(i in c(1:3)){
    gene_mhc1_filename<- paste0("temp/netMHCPan40/output/",mut_tmp,i,'.xls') 
    gene_mhc1<- read.table(gene_mhc1_filename,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
    # Only select peptides with actual mutation present
    idx_core<- which(gene_mhc1$Pos%in%0:10) # Cannever start at position 11 (12) or higher
    idx_core<- intersect(idx_core,which(nchar(gene_mhc1$Peptide)+gene_mhc1$Pos>10)) # Remove pos 0 (1) for 8/9/10-mers, ...
    gene_mhc1<- gene_mhc1[idx_core,]
    # Get ranks
    gene_mhc1_rank<- gene_mhc1[,grep("Rank",colnames(gene_mhc1))]
    colnames(gene_mhc1_rank)<- setdiff(unlist(strsplit(readLines(gene_mhc1_filename,n=1),"\t")),"")
    # Calculate PBR
    if(i==1) PBR_tmp<- apply(gene_mhc1_rank,MARGIN = 2, "min")
    else PBR_tmp<- c(PBR_tmp,apply(gene_mhc1_rank,MARGIN = 2, "min"))
  }
  
  # Put to HLA matrix
  if(m==1) PBR<- as.data.frame(matrix(NA,length(driver_muts),length(PBR_tmp),dimnames = list(driver_muts,names(PBR_tmp))))
  PBR[mut_tmp,names(PBR_tmp)]<- PBR_tmp 
}
PBR<- na.omit(PBR)
vp_mut_matrix<- vp_mut_matrix[,rownames(PBR)]

################
# Calculate ORs
################

allele_PBR_ls<- list()

for(a in colnames(PBR)){
  
  cat(which(colnames(PBR)==a),"\n")
  
  # Put allele in matrix
  vp_HLA1_mut_matrix<- vp_mut_matrix
  vp_HLA1_mut_matrix[,]<- NA
  for(i in 1:nrow(vp_HLA1_mut_matrix)){
    vp_HLA1_mut_matrix[i,]<- PBR[,a]
  }
  
  # Calculate LR models
  logreg_wp<- do_logreg(vp_mut_matrix,vp_HLA1_mut_matrix,"wp")

  # Data in list
  allele_PBR_ls[[a]]<- logreg_wp
}
saveRDS(allele_PBR_ls,"results/data/TCGA_alleles_OR.rds")

#########################################################
# Get PBRs per allele for driver muts (TCGA genotypes)
#########################################################

# Load
mut_matrix<- readRDS("data/mut_matrix.rds")
HLA_matrix<- readRDS("data/HLA_matrix.rds")
load("../immunoediting_2019/data/TCGA_HLA_types.RData")

# Get common patients
rownames(mut_matrix)<- substr(rownames(mut_matrix),1,12) 
TCGA_HLA_types_netMHC<- na.omit(TCGA_HLA_types_netMHC)
TCGA_HLA_types_netMHC<- TCGA_HLA_types_netMHC[-union(grep("HLA-B13:07",TCGA_HLA_types_netMHC[,3]),grep("HLA-B13:07",TCGA_HLA_types_netMHC[,4])),]
common_pt<- intersect(rownames(mut_matrix),rownames(TCGA_HLA_types_netMHC))
mut_matrix<- mut_matrix[common_pt,]
TCGA_HLA_types_netMHC<- TCGA_HLA_types_netMHC[common_pt,]
HLA_matrix<- HLA_matrix[common_pt,]

# Median PBRs
PBR_med_A<- median(as.numeric(unlist(PBR[grep("HLA-A",names(PBR))])),na.rm=T) # 3.42
PBR_med_B<- median(as.numeric(unlist(PBR[grep("HLA-B",names(PBR))])),na.rm=T) # 3.70
PBR_med_C<- median(as.numeric(unlist(PBR[grep("HLA-C",names(PBR))])),na.rm=T) # 2.88
PBR_med_all<- median(as.numeric(unlist(PBR)),na.rm=T) # 3/55
PHBR_med_all<- median(as.numeric(unlist(HLA_matrix)),na.rm=T) # 1.58

muts<- c(readRDS("data/driver_muts13.rds"),"PIK3CA_H1047R","TP53_R273C","IDH1_R132C")
PBR_t<- NULL
PHBR_t<- NULL
for(mut in muts){
  
  isMut<- mut_matrix[,mut]==1
  PBR_tmp<- PBR[mut,]
  PHBR_tmp<- HLA_matrix[,mut]
  for(HLA_idx in list(1:2,3:4,5:6)){
    alleles_mut<- as.character(TCGA_HLA_types_netMHC[isMut,HLA_idx])
    alleles_noMut<- as.character(TCGA_HLA_types_netMHC[!isMut,HLA_idx])
    PBR_t<- rbind(PBR_t,
                  data.frame(PBR=as.numeric(PBR_tmp[alleles_mut]),isMut=T,mut=mut,HLA=c("A","A","B","B","C","C")[HLA_idx][1]),
                  data.frame(PBR=as.numeric(PBR_tmp[alleles_noMut]),isMut=F,mut=mut,HLA=c("A","A","B","B","C","C")[HLA_idx][1])
    )
  }
  PHBR_t<- rbind(PHBR_t,
                data.frame(PHBR=as.numeric(PHBR_tmp[isMut]),isMut=T,mut=mut),
                data.frame(PHBR=as.numeric(PHBR_tmp[!isMut]),isMut=F,mut=mut)
  )
  
}

# Save PBR table
save(PBR_t,PHBR_t,PBR_med_A,PBR_med_B,PBR_med_C,PBR_med_all,PHBR_med_all,file="results/data/TCGA_alleles_PBR_t.RData")
