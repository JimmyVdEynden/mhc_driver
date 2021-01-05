# Run netMHCPan
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -1 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input_wt/","driver_mut",".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_wt/","driver_mut_TCGA","1.xls > /dev/null"))
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -2 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input_wt/","driver_mut",".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_wt/","driver_mut_TCGA","2.xls > /dev/null"))
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -3 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input_wt/","driver_mut",".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_wt/","driver_mut_TCGA","3.xls > /dev/null"))

# Functions
source("scripts/functions/harmonic_mean.R")

# Prototypical alleles
HLA_proto_mhc1<- as.character(read.table("../mhc2/data/mhc1_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))

# Get drive rmutation
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))#

# Import results
PHBR_TCGA_mhc1_wt<- rep(NA,length(driver_muts))
names(PHBR_TCGA_mhc1_wt)<- driver_muts
for(i in c(1:3)){
  gene_mhc1_filename<- paste0("temp/netMHCPan40/output_wt/driver_mut_TCGA",i,'.xls')
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
  PBR_tmp<- apply(gene_mhc1_rank,MARGIN = 2, function(x) tapply(x, gene_mhc1$ID,"min")) 
  if(i==1) PBR<- PBR_tmp
  else PBR<- cbind(PBR, PBR_tmp)
}

# PHBR for each sample
load("~/projects/immunoediting_2019/data/TCGA_HLA_types.RData")
HLA_wt_matrix<- readRDS("data/HLA_matrix.rds")
HLA_wt_matrix[,]<- NA 

for(mut_tmp in rownames(PBR)){
  cat(mut_tmp," ")
  PBR_tmp<- PBR[mut_tmp,]
  
  PBR_all_alleles<- cbind(
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-A1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-A2"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-B1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-B2"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-C1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-C2"]]
  )
  rownames(PBR_all_alleles)<- rownames(TCGA_HLA_types_netMHC)
  PHBR_tmp<- apply(PBR_all_alleles,1,"harmonic_mean")
  
  # Put to HLA matrix
  common_samples<- intersect(rownames(HLA_wt_matrix),names(PHBR_tmp)) # 8760
  HLA_wt_matrix[common_samples,mut_tmp]<- PHBR_tmp[common_samples]
  
}

saveRDS(HLA_wt_matrix,file="data/HLA_wt_matrix.rds")
