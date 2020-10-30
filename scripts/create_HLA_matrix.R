# Load mut matrix
mut_matrix<- readRDS("data/mut_matrix.rds")

# Load driver mutation data (needed to get ENSP id)
load("temp/TCGA_maf_driver.RData")

# Harmonic mean function
source("scripts/functions/harmonic_mean.R")

# Load TCGA Polysolver HLA genotypes & save all alleles to txt file to use with netMHCPan (max 1024 chars per line!)
load("~/projects/immunoediting_2019/data/TCGA_HLA_types.RData")
HLA_alleles_TCGA<- sort(setdiff(na.omit(unique(c(TCGA_HLA_types_netMHC[,1],TCGA_HLA_types_netMHC[,2],TCGA_HLA_types_netMHC[,3],TCGA_HLA_types_netMHC[,4],TCGA_HLA_types_netMHC[,5],TCGA_HLA_types_netMHC[,6]))),""))
HLA_alleles_TCGA<- setdiff(HLA_alleles_TCGA,"HLA-B13:07") # B13:07 not in netMHCPan!
cat(HLA_alleles_TCGA[1:65],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt")
cat("\n",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
cat(HLA_alleles_TCGA[66:130],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
cat("\n",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
cat(HLA_alleles_TCGA[131:length(HLA_alleles_TCGA)],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)

# Create empty HLA matrix
HLA_matrix<- mut_matrix
HLA_matrix[,]<- NA

# Get peptides for driver muts
library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information
edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) # Restrict edb to coding genes
for(mut_tmp in colnames(HLA_matrix)[1:ncol(HLA_matrix)]){
  
  cat(grep(mut_tmp,colnames(HLA_matrix)),"\n")
  mut_tmp_gene<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Hugo_Symbol"])
  mut_tmp_ENSP<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"ENSP"])
  mut_tmp_aa<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Amino_acids"])
  mut_tmp_aa_ref<- gsub("\\/.","",mut_tmp_aa)
  mut_tmp_aa_alt<- gsub(".\\/","",mut_tmp_aa)
  mut_tmp_pos<- as.numeric(unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Protein_position"]))
  
  # Get protein sequence information
  if(mut_tmp_gene=="RUNX1T1") mut_tmp_ENSP<- "ENSP00000428742" # Manual correction, mutations only matched this ENSP
  if(mut_tmp_gene=="CRLF2") mut_tmp_ENSP<- "ENSP00000370978" # Manual correction, mutations only matched this ENSP
  if(mut_tmp_gene=="TCF7L2") mut_tmp_ENSP<- "ENSP00000486891" # Manual correction, mutations only matched this ENSP
  if(mut_tmp_gene=="MYC") mut_tmp_ENSP<- "ENSP00000478887" # Manual correction, mutations only matched this ENSP
  if(mut_tmp_ENSP=="ENSP00000298047") next # not found
  prot_sequence<- proteins(edb, filter = ProteinIdFilter(mut_tmp_ENSP), return.type = "AAStringSet")
  # Use gene if ENSP not found
  if(length(prot_sequence)==0){ 
    proteins_for_gene<- proteins(edb, filter = GeneNameFilter(mut_tmp_gene), return.type = "AAStringSet")
    prot_sequence <- proteins_for_gene[which(width(proteins_for_gene)==max(width(proteins_for_gene)))[1]] # Take max if multiple
  }
  
  # Get 21-mers
  sequence_strings <- as.character(prot_sequence)
  pep11_tmp<- substr(sequence_strings,mut_tmp_pos-10,mut_tmp_pos+10)
  # IGV check: ok for BRAF V600E
  
  # Add mutations
  pos_in_pep_tmp<- 11
  if(mut_tmp_pos<11) pos_in_pep_tmp<-  mut_tmp_pos
  if(substr(pep11_tmp,pos_in_pep_tmp,pos_in_pep_tmp)!=mut_tmp_aa_ref) stop("Check mutation: ", mut_tmp)
  # Adapt the following after errors: "APCS2307L"
  substr(pep11_tmp,pos_in_pep_tmp,pos_in_pep_tmp)<- mut_tmp_aa_alt
  
  # Save peptide faste files
  cat(">",mut_tmp,"\n",pep11_tmp,"\n",sep="",file = paste0("temp/netMHCPan40/input/",mut_tmp,".fasta"))
  
  # Run NetMHCPan: split in 3 because netMHCPan only accepts HLA allele file with fewer than 1024 chars
  system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -1 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"1.xls > /dev/null"))
  system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -2 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"2.xls > /dev/null"))
  system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -3 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"3.xls > /dev/null"))

  # Import results
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
  
  # PHBR for each sample
  PBR_all_alleles<- cbind(
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-A1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-A2"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-B1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-B2"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-C1"]],
    PBR_tmp[TCGA_HLA_types_netMHC[,"HLA-C2"]]
  )
  rownames(PBR_all_alleles)<- rownames(TCGA_HLA_types_netMHC)
  PHBR_tmp<- apply(PBR_all_alleles,1,"harmonic_mean") # 9.566061 for prototypical (1.91 for wt)
  
  # Put to HLA matrix
  rownames(HLA_matrix)<- substr(rownames(HLA_matrix),1,12)
  common_samples<- intersect(rownames(HLA_matrix),names(PHBR_tmp)) # 8760
  HLA_matrix[common_samples,mut_tmp]<- PHBR_tmp[common_samples]
}

# Save
saveRDS(HLA_matrix,file="data/HLA_matrix.rds")

