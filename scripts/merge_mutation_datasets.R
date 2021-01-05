############################
# Merge mutation datasets
############################

# Get genome coordinates & PHBR
for(mut_type in c("","_recurr_noDriver","_noRecurr_driver","_noRecurr_noDriver")){
  
  cat(mut_type,"\n")
  isDriver<- !grepl("noDriver",mut_type)
  isRecurr<- !grepl("noRecurr",mut_type)
  
  # mut
  if(mut_type=="") load("data/TCGA_maf_driver.RData")
  else load(paste0("data/TCGA_maf",mut_type,".RData"))
  
  # Only unique
  TCGA_maf_driver<- TCGA_maf_driver[!duplicated(TCGA_maf_driver$mut_id),]
  rownames(TCGA_maf_driver)<- TCGA_maf_driver$mut_id
  
  # Get genome coordinates
  TCGA_maf_driver<- TCGA_maf_driver[,1:4]
  
  # Add prototype PHBR
  TCGA_maf_driver$PHBR<- NA
  library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information
  edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) 
  
  for(mut_tmp in rownames(TCGA_maf_driver)){
    #PHBR
    for(i in c(1:3)){
      gene_mhc1_filename<- paste0("temp/netMHCPan40/output",mut_type,"/",mut_tmp,i,'.xls') 
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
    
    # PHBR for each sample
    source("scripts/functions/harmonic_mean.R")
    HLA_proto<- as.character(read.table("../mhc2/data/mhc1_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))
    PBR_all_alleles<- PBR_tmp[HLA_proto]   
    TCGA_maf_driver[mut_tmp,"PHBR"]<- harmonic_mean(PBR_all_alleles) # 9.566061 for prototypical (1.91 for wt)
  } 
  
  # Fuse with other data
  if(mut_type=="") driver_overview<- cbind(TCGA_maf_driver,isCGC=isDriver,isRecurr=isRecurr)
  else driver_overview<- rbind(driver_overview,cbind(TCGA_maf_driver,isCGC=isDriver,isRecurr=isRecurr))
}

# Add data for unobserved mutations
####################################

# 15000 random mutations (hg38!)
GPPM_rand <- readRDS("data/GPPM_rand.rds")

# TCGA mutations
TCGA_maf<- readRDS("data/TCGA_maf_MC3_selection.rds")
TCGA_maf$mut_id<- paste0(TCGA_maf$Hugo_Symbol,"_",gsub("p\\.","",TCGA_maf$HGVSp_Short))

# Get 500 random mutations that are never observed
mut_zero<- setdiff(GPPM_rand$mut,TCGA_maf$mut_id) # 14589
idx_CGC<- which(gsub("\\_.*","",mut_zero)%in%readRDS("downloads/cosmic/CGC/v91/CGC_v91.rds"))
mut_zero<- mut_zero[-idx_CGC] # No CGC genes
GPPM_rand_sel<- GPPM_rand[GPPM_rand$mut%in%mut_zero]
GPPM_rand_sel<- GPPM_rand_sel[!is.na(GPPM_rand_sel$PHBR_mut),] # PHBR not NA
GPPM_rand_sel<- GPPM_rand_sel[sample(1:length(GPPM_rand_sel),500)] # 500 random muts

# Data frame
mut_zero_df<- data.frame(
  Hugo_Symbol=GPPM_rand_sel$gene,
  Chromosome=seqnames(GPPM_rand_sel),
  Start_Position=pos(GPPM_rand_sel),
  End_Position=pos(GPPM_rand_sel),
  PHBR=GPPM_rand_sel$PHBR_mut,
  isCGC=F,
  isRecurr=F,
  isObserved=F,
  row.names = mut_zero_df$mut)

# liftOver hg38 to hg19
# cat(paste0("chr",mut_zero_df$Chromosome,":",mut_zero_df$Start_Position,"-",mut_zero_df$End_Position),sep="\n")
# manually in https://genome.ucsc.edu/cgi-bin/hgLiftOver and upload file
lo_res<- read.table("temp/hglft_genome_3b0d5_b9ca60.bed")
lo_res_pos<- as.numeric(gsub(".*\\-","",lo_res$V1))
mut_zero_df$Start_Position<- lo_res_pos
mut_zero_df$End_Position<- lo_res_pos

# Merge with overview
driver_overview$isObserved<- T
driver_overview<- rbind(driver_overview,mut_zero_df)

######################################
# Add evolutionary conservation scores
######################################

library(GenomicScores)
library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19
phastCon<- gscores(phast, GRanges(seqnames=paste0("chr",driver_overview$Chromosome), IRanges(start=driver_overview$Start_Position-10, end=driver_overview$Start_Position+10)))
driver_overview$phastCons100way<- phastCon$default

phylo<- getGScores("phyloP100way.UCSC.hg19") # Takes some time first time, then faster
# gscores(phylo, GRanges(seqnames="chr22", IRanges(start=50967020:50967025, width=1)))
# phyloP<- gscores(phylo, GRanges(seqnames=paste0("chr",driver_overview$Chromosome), IRanges(start=driver_overview$Start_Position, width=1)))
phyloP<- gscores(phylo, GRanges(seqnames=paste0("chr",driver_overview$Chromosome), IRanges(start=driver_overview$Start_Position-10, width=21)))
driver_overview$phyloP100way<- phyloP$default

# Save
saveRDS(driver_overview,file="data/driver_overview_cons.rds")
