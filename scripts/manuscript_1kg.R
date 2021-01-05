source("scripts/functions/harmonic_mean.R")
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")
muts<- readRDS("data/driver_muts13.rds")
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))
source("scripts/functions/do_logreg.R")

# Load genotype data 1kg
pop_HLA<- read.table("downloads/1kg/1000_genomes_hla.tsv",sep="\t",header = T,row.names = 3,stringsAsFactors = F)
idx_amb<- unique(as.numeric(unlist(apply(pop_HLA[,3:8],2,function(x) grep("\\/|\\*|None",x))))) # identify ambiguous samples (e.g. 02/03 format - * labelled - "none"): 116  in total
pop_HLA<- pop_HLA[-idx_amb,] 
pops<- pop_HLA$Region
pop_HLA<- pop_HLA[,3:8] # Focus on MHC-I only (MHC-II not complete)

table(pops)
# AFR AMR EAS EUR SAS 
# 674 359 511 512 521 

####################################
# (1) PROTOTYPES PER POPULATION
####################################

# Get prototypes per population
##################################
proto_1kg<- as.data.frame(matrix(NA,length(unique(pops)),6,dimnames = list(unique(pops),c("A1","A2","B1","B2","C1","C2"))))
proto_1kg_freq<- proto_1kg
for(pop in unique(pops)){
  pop_HLA_sel<- pop_HLA[pops==pop,] 
  freq_A<- sort(prop.table(table(c(pop_HLA_sel$HLA.A.1,pop_HLA_sel$HLA.A.2))),decreasing = T)[1:2]
  names(freq_A)<- paste0("HLA-A",names(freq_A))
  freq_B<- sort(prop.table(table(c(pop_HLA_sel$HLA.B.1,pop_HLA_sel$HLA.B.2))),decreasing = T)[1:2]
  names(freq_B)<- paste0("HLA-B",names(freq_B))
  freq_C<- sort(prop.table(table(c(pop_HLA_sel$HLA.C.1,pop_HLA_sel$HLA.C.2))),decreasing = T)[1:2]
  names(freq_C)<- paste0("HLA-C",names(freq_C))
  proto_1kg[pop,]<- c(names(freq_A),names(freq_B),names(freq_C))
  proto_1kg_freq[pop,]<- c(freq_A,freq_B,freq_C)
}
proto_1kg
# A1         A2         B1         B2         C1         C2
# AFR HLA-A02:01 HLA-A23:01 HLA-B53:01 HLA-B35:01 HLA-C04:01 HLA-C16:01
# AMR HLA-A02:01 HLA-A24:02 HLA-B44:03 HLA-B07:02 HLA-C04:01 HLA-C07:02
# EAS HLA-A11:01 HLA-A24:02 HLA-B46:01 HLA-B40:01 HLA-C01:02 HLA-C07:02
# EUR HLA-A02:01 HLA-A03:01 HLA-B07:02 HLA-B44:02 HLA-C04:01 HLA-C07:01
# SAS HLA-A11:01 HLA-A01:01 HLA-B40:06 HLA-B52:01 HLA-C06:02 HLA-C07:02
proto_1kg_freq
# AFR 0.1142433 0.09347181 0.16246291 0.08976261 0.2462908 0.1201780
# AMR 0.2200557 0.14066852 0.07381616 0.05849582 0.1922006 0.1100279
# EAS 0.2338552 0.19373777 0.15264188 0.10958904 0.1937378 0.1712329
# EUR 0.3007812 0.15625000 0.11621094 0.08691406 0.1318359 0.1250000
# SAS 0.1660269 0.16314779 0.09596929 0.09213052 0.1420345 0.1238004

# Calculate matrices & log reg models
######################################

pop_PHBR_ls<- list()

for(pop in rownames(proto_1kg)){
  
  cat(pop," ...","\n")
  # Get proto HLA genotype
  HLA_proto_mhc1<- as.character(proto_1kg[pop,])
  
  # Calculate PHBR
  PHBR_proto_mhc1_mut<- rep(NA,length(driver_muts))
  names(PHBR_proto_mhc1_mut)<- driver_muts
  for(mut_tmp in driver_muts){
    for(i in c(1:3)){
      gene_mhc1_filename<- paste0("temp/netMHCPan40/output/",mut_tmp,i,'.xls')
      if(!file.exists(gene_mhc1_filename)) next
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
    if(!file.exists(gene_mhc1_filename)) next
    
    PBR_tmp<- PBR_tmp[HLA_proto_mhc1]
    PHBR_proto_mhc1_mut[mut_tmp]<- harmonic_mean(PBR_tmp) # 9.566061 for prototypical (1.91 for wt)
  }
  
  # Put proto in matrix
  vp_HLA1_mut_matrix<- vp_mut_matrix
  vp_HLA1_mut_matrix[,]<- NA
  common_muts<- intersect(names(PHBR_proto_mhc1_mut),colnames(vp_HLA1_mut_matrix))
  for(i in 1:nrow(vp_HLA1_mut_matrix)){
    vp_HLA1_mut_matrix[i,common_muts]<- PHBR_proto_mhc1_mut[common_muts]
  }
  
  # Calculate LR models
  logreg_wp<- do_logreg(vp_mut_matrix,vp_HLA1_mut_matrix,"wp")
  idx_mut13<- which(colnames(vp_mut_matrix)%in%muts)
  logreg_wp_excl<- do_logreg(vp_mut_matrix[,-idx_mut13],vp_HLA1_mut_matrix[,-idx_mut13],"wp")
  logreg<- rbind(wp=logreg_wp,wp_excl=logreg_wp_excl)
  
  # Data in list
  pop_PHBR_ls[[pop]]<- list(vp_HLA1_mut_matrix=vp_HLA1_mut_matrix,logreg=logreg)
}
saveRDS(pop_PHBR_ls,"results/data/1kg_LR.rds")
saveRDS(proto_1kg,"results/data/1kg_proto.rds")

#############################################
# (2) LOWEST AF GENOTYPES PER POPULATION
#############################################

# Get genotypes for lowest AF per population
##############################################

# Get TCGA alleles
TCGA_alleles<- NULL
for(i in 1:3){
  gene_mhc1_filename<- paste0("temp/netMHCPan40/output/BRAF_V600E",i,'.xls')
  TCGA_alleles<- c(TCGA_alleles,setdiff(unlist(strsplit(readLines(gene_mhc1_filename,n=1),"\t")),""))
}

# Get random alleles that were never called in a population
proto_1kg<- as.data.frame(matrix(NA,length(unique(pops)),6,dimnames = list(unique(pops),c("A1","A2","B1","B2","C1","C2"))))
for(pop in unique(pops)){
  pop_HLA_sel<- pop_HLA[pops==pop,] 
  HLA_A<- sample(setdiff(TCGA_alleles[grep("HLA-A",TCGA_alleles)],paste0("HLA-A",c(pop_HLA_sel$HLA.A.1,pop_HLA_sel$HLA.A.2))),2)
  HLA_B<- sample(setdiff(TCGA_alleles[grep("HLA-B",TCGA_alleles)],paste0("HLA-B",c(pop_HLA_sel$HLA.B.1,pop_HLA_sel$HLA.B.2))),2)
  HLA_C<- sample(setdiff(TCGA_alleles[grep("HLA-C",TCGA_alleles)],paste0("HLA-C",c(pop_HLA_sel$HLA.C.1,pop_HLA_sel$HLA.C.2))),2)
  proto_1kg[pop,]<- c(HLA_A,HLA_B,HLA_C)
}
proto_1kg
#         A1         A2         B1         B2         C1         C2
# AFR HLA-A11:02 HLA-A02:07 HLA-B15:20 HLA-B55:04 HLA-C07:05 HLA-C16:02
# AMR HLA-A02:03 HLA-A69:01 HLA-B15:58 HLA-B15:12 HLA-C07:17 HLA-C18:01
# EAS HLA-A11:04 HLA-A33:01 HLA-B15:08 HLA-B39:10 HLA-C18:01 HLA-C08:04
# EUR HLA-A11:02 HLA-A02:10 HLA-B15:21 HLA-B73:01 HLA-C14:03 HLA-C07:05
# SAS HLA-A30:04 HLA-A29:02 HLA-B27:03 HLA-B73:01 HLA-C04:04 HLA-C02:10


# Calculate matrices & log reg models
######################################

pop_PHBR_ls<- list()

for(pop in rownames(proto_1kg)){
  
  cat(pop," ...","\n")
  # Get proto HLA genotype
  HLA_proto_mhc1<- as.character(proto_1kg[pop,])
  
  # Calculate PHBR
  PHBR_proto_mhc1_mut<- rep(NA,length(driver_muts))
  names(PHBR_proto_mhc1_mut)<- driver_muts
  for(mut_tmp in driver_muts){
    for(i in c(1:3)){
      gene_mhc1_filename<- paste0("temp/netMHCPan40/output/",mut_tmp,i,'.xls')
      if(!file.exists(gene_mhc1_filename)) next
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
    if(!file.exists(gene_mhc1_filename)) next
    
    PBR_tmp<- PBR_tmp[HLA_proto_mhc1]
    PHBR_proto_mhc1_mut[mut_tmp]<- harmonic_mean(PBR_tmp) # 9.566061 for prototypical (1.91 for wt)
  }
  
  # Put proto in matrix
  vp_HLA1_mut_matrix<- vp_mut_matrix
  vp_HLA1_mut_matrix[,]<- NA
  common_muts<- intersect(names(PHBR_proto_mhc1_mut),colnames(vp_HLA1_mut_matrix))
  for(i in 1:nrow(vp_HLA1_mut_matrix)){
    vp_HLA1_mut_matrix[i,common_muts]<- PHBR_proto_mhc1_mut[common_muts]
  }
  
  # Calculate LR models
  logreg_wp<- do_logreg(vp_mut_matrix,vp_HLA1_mut_matrix,"wp")
  idx_mut13<- which(colnames(vp_mut_matrix)%in%muts)
  logreg_wp_excl<- do_logreg(vp_mut_matrix[,-idx_mut13],vp_HLA1_mut_matrix[,-idx_mut13],"wp")
  logreg<- rbind(wp=logreg_wp,wp_excl=logreg_wp_excl)
  
  # Data in list
  pop_PHBR_ls[[pop]]<- list(vp_HLA1_mut_matrix=vp_HLA1_mut_matrix,logreg=logreg)
}
saveRDS(pop_PHBR_ls,"results/data/1kg_LR_lowFreq.rds")
saveRDS(proto_1kg,"results/data/1kg_proto_lowFreq.rds")

