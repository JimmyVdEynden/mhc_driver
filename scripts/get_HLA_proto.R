# Functions
source("scripts/functions/harmonic_mean.R")

# Prototypical alleles
HLA_proto_mhc1<- as.character(read.table("../mhc2/data/mhc1_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))
HLA_proto_mhc2<- as.character(read.table("../mhc2/data/mhc2_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))

# Get driver mutations
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))

# Create MHC-I prototypical PHBR
################################

# CALCULATE MHC1 WT AFFINITIES FOR PROTO

# Load driver mutation data (needed to get ENSP id)
load("temp/TCGA_maf_driver.RData")

# Get peptides for driver muts
library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information: add to yml script
edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) # Restrict edb to coding genes
for(mut_tmp in driver_muts){

  i<- grep(mut_tmp,driver_muts)
  cat(i,"\n")
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
  prot_sequence <- proteins(edb, filter = ProteinIdFilter(mut_tmp_ENSP), return.type = "AAStringSet")
  if(length(prot_sequence)==0){
    proteins_for_gene <- proteins(edb, filter = GeneNameFilter(mut_tmp_gene), return.type = "AAStringSet")
    prot_sequence<- proteins_for_gene[which(width(proteins_for_gene)==max(width(proteins_for_gene)))[1]] # Take longest isoform if ENSP not present
  }

  sequence_strings <- as.character(prot_sequence)
  pep11_tmp<- substr(sequence_strings,mut_tmp_pos-10,mut_tmp_pos+10)
  # IGV check: ok for BRAF V600E

  # Save peptide faste files
  if(i==1) cat(">",mut_tmp,"\n",pep11_tmp,"\n",sep="",file = paste0("temp/netMHCPan40/input_wt/","driver_mut",".fasta"))
  else cat(">",mut_tmp,"\n",pep11_tmp,"\n",sep="",file = paste0("temp/netMHCPan40/input_wt/","driver_mut",".fasta"),append=T)
}

# Run NetMHCPan
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head ../mhc2/data/mhc1_HLA_alleles_proto.txt` -f temp/netMHCPan40/input_wt/","driver_mut",".fasta -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_wt/","driver_mut",".xls > /dev/null"))

# Import results
gene_mhc1_filename_wt<- paste0("temp/netMHCPan40/output_wt/driver_mut.xls")
if(!file.exists(gene_mhc1_filename_wt)) next
gene_mhc1_wt<- read.table(gene_mhc1_filename_wt,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
# Only select peptides with actual mutation present
idx_core<- which(gene_mhc1_wt$Pos%in%0:10) # Can never start at position 11 (12) or higher
idx_core<- intersect(idx_core,which(nchar(gene_mhc1_wt$Peptide)+gene_mhc1_wt$Pos>10)) # Remove pos 0 (1) for 8/9/10-mers, ...
gene_mhc1_wt<- gene_mhc1_wt[idx_core,]
# Get ranks
gene_mhc1_wt_rank<- gene_mhc1_wt[,grep("Rank",colnames(gene_mhc1_wt))]
colnames(gene_mhc1_wt_rank)<- setdiff(unlist(strsplit(readLines(gene_mhc1_filename_wt,n=1),"\t")),"")

# Calculate PBR & PHBR
PBR_tmp<- apply(gene_mhc1_wt_rank,MARGIN = 2, function(x) tapply(x, gene_mhc1_wt$ID,"min"))
PHBR_proto_mhc1_wt<- apply(PBR_tmp,1,harmonic_mean)
PHBR_proto_mhc1_wt<- PHBR_proto_mhc1_wt[driver_muts]

# MHC 1 MUT AFFINITIES

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

# Create MHC-II prototypical PHBR
##################################
PHBR_proto_mhc2_mut<- rep(NA,length(driver_muts))
names(PHBR_proto_mhc2_mut)<- driver_muts

# CALCULATE MHC2 MUT AFFINITIES FOR PROTO

# Load driver mutation data (needed to get ENSP id)
load("temp/TCGA_maf_driver.RData")

# Get peptides for driver muts
library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information: add to yml script
edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) # Restrict edb to coding genes
for(mut_tmp in driver_muts){

  i<- grep(mut_tmp,driver_muts)
  cat(i,"\n")
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
  prot_sequence <- proteins(edb, filter = ProteinIdFilter(mut_tmp_ENSP), return.type = "AAStringSet")
  if(length(prot_sequence)==0){
    proteins_for_gene <- proteins(edb, filter = GeneNameFilter(mut_tmp_gene), return.type = "AAStringSet")
    prot_sequence<- proteins_for_gene[which(width(proteins_for_gene)==max(width(proteins_for_gene)))[1]] # Take longest isoform if ENSP not present
  }

  sequence_strings <- as.character(prot_sequence)
  pep15_tmp<- substr(sequence_strings,mut_tmp_pos-14,mut_tmp_pos+14)
  # IGV check: ok for BRAF V600E

  # Add mutations
  pos_in_pep_tmp<- 15
  if(mut_tmp_pos<pos_in_pep_tmp) pos_in_pep_tmp<-  mut_tmp_pos
  if(substr(pep15_tmp,pos_in_pep_tmp,pos_in_pep_tmp)!=mut_tmp_aa_ref) stop("Check mutation: ", mut_tmp)
  substr(pep15_tmp,pos_in_pep_tmp,pos_in_pep_tmp)<- mut_tmp_aa_alt

  # Save peptide faste files
  if(i==1) cat(">",mut_tmp,"\n",pep15_tmp,"\n",sep="",file = paste0("temp/netmhciipan32/input/","driver_mut",".fasta"))
  else cat(">",mut_tmp,"\n",pep15_tmp,"\n",sep="",file = paste0("temp/netmhciipan32/input/","driver_mut",".fasta"),append=T)
}

# Run NetMHCPan
system(paste0("/home/labgroups/ccgg/tools/netMHCIIpan-3.2/netMHCIIpan -a `head /home/jimmy/projects/mhc2/data/mhc2_HLA_alleles_proto.txt` -f temp/netmhciipan32/input/","driver_mut",".fasta -inptype 0 -length 15 -xls -xlsfile temp/netmhciipan32/output/","driver_mut",".xls > /dev/null"))

# Import results
gene_mhc2_filename<- paste0("temp/netmhciipan32/output/","driver_mut",'.xls') 
gene_mhc2<- read.table(gene_mhc2_filename,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
gene_mhc2_rank<- gene_mhc2[,grep("Rank",colnames(gene_mhc2))]
colnames(gene_mhc2_rank)<- setdiff(unlist(strsplit(readLines(gene_mhc2_filename,n=1),"\t")),"")

# Calculate PBR & PHBR
PBR_proto<- apply(gene_mhc2_rank,MARGIN = 2, function(x) tapply(x, gene_mhc2$ID,"min"))
PHBR_proto_mhc2_mut<- apply(PBR_proto,1,harmonic_mean)
PHBR_proto_mhc2_mut<- PHBR_proto_mhc2_mut[driver_muts]
  
# Save
######
save(PHBR_proto_mhc1_mut,PHBR_proto_mhc1_wt,PHBR_proto_mhc2_mut,file="data/PHBR_proto.RData")

