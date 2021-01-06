###########################
# Mut matrix
###########################

# Mutations
CCLE_mutations<- readr::read_csv("downloads/CCLE/CCLE_mutations.csv")

# Only missense
CCLE_maf<- CCLE_mutations[CCLE_mutations$Variant_Classification=="Missense_Mutation",]

# Factorize samples
CCLE_maf$DepMap_ID<- factor(CCLE_maf$DepMap_ID) # 1749 samples

# Add mut_id
aa_change<- gsub("p\\.","",CCLE_maf$Protein_Change)
CCLE_maf$mut_id<- paste0(CCLE_maf$Hugo_Symbol,"_",aa_change)

# Get driver muts from other study
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))

# Select these mutation data: 1500 muts
CCLE_maf_driver<- subset(CCLE_maf,CCLE_maf$mut_id%in%driver_muts) 
CCLE_maf_driver$mut_id<- factor(CCLE_maf_driver$mut_id,levels = driver_muts)

# Create matrix
mut_matrix<- table(CCLE_maf_driver$DepMap_ID,CCLE_maf_driver$mut_id)
# dim(mut_matrix) # 1749 x 688

# Sort driver muts by frequency
mut_matrix<- mut_matrix[,order(colSums(mut_matrix),decreasing = T)]
colSums(mut_matrix)[1:10]
# BRAF_V600E     KRAS_G12D    TP53_R273H     KRAS_G12V    TP53_R175H    TP53_R273C     NRAS_Q61K PIK3CA_H1047R 
# 112            71            48            46            44            39            38            33 
# TP53_R248W     KRAS_G12C 
# 31            28 


###########################
# HLA matrix
###########################

# Sample info
sample_info<- as.data.frame(readr::read_csv("downloads/CCLE/sample_info.csv"))
sample_info<- sample_info[!duplicated(sample_info$stripped_cell_line_name)&!is.na(sample_info$stripped_cell_line_name),]
rownames(sample_info)<- sample_info$stripped_cell_line_name

# HLA data From TRON project
HLA<- as.data.frame(readr::read_csv("downloads/CCLE/TRON_seq2hla.csv"))
rownames(HLA)<- toupper(gsub("\\.|\\,|\\-|\\_","",HLA$CELLLINE))
rownames(HLA)[rownames(HLA)=="D341MED"]<- "D341" 

# Get HLA alleles
HLA_t<- data.frame(
  HLA_A1=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-A-A1`)),
  HLA_A2=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-A-A2`)),
  HLA_B1=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-B-A1`)),
  HLA_B2=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-B-A2`)),
  HLA_C1=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-C-A1`)),
  HLA_C2=paste0("HLA-",gsub("\\*|\\'","",HLA$`HLA-C-A2`)),stringsAsFactors = F, row.names = sample_info[rownames(HLA),"DepMap_ID"]
)


# Generate HLA matrix
########################

# Harmonic mean function
source("scripts/functions/harmonic_mean.R")

# Create empty HLA matrix
driver_muts<- colnames(mut_matrix)
HLA_matrix<- mut_matrix
HLA_matrix[,]<- NA

# HLA alleles not called bedore
HLA_alleles_called<- c(
  setdiff(unlist(strsplit(readLines(paste0("temp/netMHCPan40/output/BRAF_V600E",1,'.xls'),n=1),"\t")),""),
  setdiff(unlist(strsplit(readLines(paste0("temp/netMHCPan40/output/BRAF_V600E",2,'.xls'),n=1),"\t")),""),
  setdiff(unlist(strsplit(readLines(paste0("temp/netMHCPan40/output/BRAF_V600E",3,'.xls'),n=1),"\t")),"")
)

HLA_all<- unique(c(
  as.character(HLA_t$HLA_A1),
  as.character(HLA_t$HLA_A2),
  as.character(HLA_t$HLA_B1),
  as.character(HLA_t$HLA_B2),
  as.character(HLA_t$HLA_C1),
  as.character(HLA_t$HLA_C2)
))

HLA_alleles_uncalled<- setdiff(HLA_all,HLA_alleles_called) # 75
HLA_alleles_uncalled<- setdiff(HLA_alleles_uncalled,"HLA-no")
cat(HLA_alleles_uncalled,sep=",",file="temp/netMHCPan40/CCLE_other_HLA_alleles.txt")

# Get PHBR scores
for(mut_tmp in driver_muts){
  
  cat(which(driver_muts==mut_tmp)," ")
  # Run NetMHCPan: for HLA alleles not called before
  if(!file.exists(paste0("temp/netMHCPan40/input/",mut_tmp,".fasta"))) next
  if(!file.exists(paste0("temp/netMHCPan40/output/",mut_tmp,"5.xls"))) system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head temp/netMHCPan40/Other_HLA_alleles.txt` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"5.xls > /dev/null"))
  
  # Import results
  for(i in c(1:5)){
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
  
  # PHBR for each sample
  PBR_all_alleles<- cbind(
    PBR_tmp[HLA_t$HLA_A1],
    PBR_tmp[HLA_t$HLA_A2],
    PBR_tmp[HLA_t$HLA_B1],
    PBR_tmp[HLA_t$HLA_B2],
    PBR_tmp[HLA_t$HLA_C1],
    PBR_tmp[HLA_t$HLA_C2]
  )
  rownames(PBR_all_alleles)<- rownames(HLA_t)
  PHBR_tmp<- apply(PBR_all_alleles,1,"harmonic_mean") 
  
  # Put to HLA matrix
  common_samples<- intersect(rownames(HLA_matrix),names(PHBR_tmp)) # 8760
  HLA_matrix[common_samples,mut_tmp]<- PHBR_tmp[common_samples]
}

save(mut_matrix, HLA_matrix, CCLE_maf_driver, file="data/CCLE_data.RData")

