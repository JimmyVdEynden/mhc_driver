# Get 100 values for each subst
################################

# Load GPPM
# GPPM<- readRDS("../mhc2/data/GPPMhm_all.rds")
# object.size(GPPM) # 88,169,300,336 bytes

# Only missense
GPPM<- GPPM[GPPM$variant=="nonsynonymous SNV"]

# selected subst
GPPM$subst<- paste0(GPPM$ref_aa,">",GPPM$alt_aa) 
subst<- unique(GPPM$subst) 
length(subst) # 150

# Select 100 random subst for each type
for(i in 1:length(subst)){
  cat(i," ")
  s<- subst[i]
  idx_sel<- sample(which(GPPM$subst==s),100)
  GPPM_tmp<- GPPM[idx_sel]
  if(i==1) GPPM_rand<- GPPM_tmp
  else GPPM_rand<- c(GPPM_rand,GPPM_tmp)
}

# Get aa position in ENSP
library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information
edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) # Restrict edb to coding genes

GPPM_rand$aa_pos<- NA
for(i in 1:length(GPPM_rand)){
  # for(i in 1:10){
  cat(i, " ")
  aa_pos_tmp<- genomeToProtein(GPPM_rand[i],edb)
  aa_pos_tmp<- start(unlist(aa_pos_tmp)[GPPM_rand[i]$ENSP])
  GPPM_rand[i]$aa_pos<- aa_pos_tmp
}

# Add 21-mer information
ENSP<- GPPM_rand$ENSP
prot_sequence<- proteins(edb, filter = ProteinIdFilter(ENSP), return.type = "AAStringSet")
sequence_strings<- as.character(prot_sequence[ENSP])
pep21<- substr(sequence_strings,GPPM_rand$aa_pos-10,GPPM_rand$aa_pos+10)

# Add mut
pep21_mut<- pep21
idx_sel<- which(GPPM_rand$aa_pos<11)
substr(pep21_mut[idx_sel],GPPM_rand[idx_sel]$aa_pos,GPPM_rand[idx_sel]$aa_pos)<- GPPM_rand[idx_sel]$alt_aa 
substr(pep21_mut[-idx_sel],11,11)<- GPPM_rand[-idx_sel]$alt_aa 

# Add pep to GPPM_rand
GPPM_rand$pep21<- pep21
GPPM_rand$pep21_mut<- pep21_mut

# Add mut to GPPM_rand
GPPM_rand$mut<- paste0(GPPM_rand$gene,"_",GPPM_rand$ref_aa,GPPM_rand$aa_pos,GPPM_rand$alt_aa)

# Check with subst
sum(substr(sequence_strings,GPPM_rand$aa_pos,GPPM_rand$aa_pos)!=GPPM_rand$ref_aa) # 0 ok

# Save
saveRDS(GPPM_rand,file="temp/GPPM_rand.rds")

# Save peptide faste files
cat(paste0(">",GPPM_rand$mut,"\n",GPPM_rand$pep21),sep="\n",file = paste0("temp/netMHCPan40/input_GPPM_rand/GPPM_rand_wt.fasta"))
cat(paste0(">",GPPM_rand$mut,"\n",GPPM_rand$pep21_mut),sep="\n",file = paste0("temp/netMHCPan40/input_GPPM_rand/GPPM_rand_mut.fasta"))

# Run NetMHCPan: split in 3 because netMHCPan only accepts HLA allele file with fewer than 1024 chars
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head ../mhc2/data/mhc1_HLA_alleles_proto.txt` -f temp/netMHCPan40/input_GPPM_rand/GPPM_rand_wt.fasta -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_GPPM_rand/GPPM_rand_wt.xls > /dev/null"))
system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head ../mhc2/data/mhc1_HLA_alleles_proto.txt` -f temp/netMHCPan40/input_GPPM_rand/GPPM_rand_mut.fasta -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output_GPPM_rand/GPPM_rand_mut.xls > /dev/null"))

# Import results wt
source("scripts/functions/harmonic_mean.R")
gene_mhc1_filename_wt<- paste0("temp/netMHCPan40/output_GPPM_rand/GPPM_rand_wt.xls")
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
# Sometimes peptides in beginning of protein and no 38 used for analysis, some by coincidence used multiple times; filter on only 38
pep_excl<- names(table(gene_mhc1_wt$ID)[table(gene_mhc1_wt$ID)!=38]) # 516
PHBR_proto_mhc1_wt[pep_excl]<- NA
# Add
GPPM_rand$PHBR_wt<- PHBR_proto_mhc1_wt[GPPM_rand$mut]

# Import results mut
gene_mhc1_filename_mut<- paste0("temp/netMHCPan40/output_GPPM_rand/GPPM_rand_mut.xls")
if(!file.exists(gene_mhc1_filename_mut)) next
gene_mhc1_mut<- read.table(gene_mhc1_filename_mut,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
# Only select peptides with actual mutation present
idx_core<- which(gene_mhc1_mut$Pos%in%0:10) # Can never start at position 11 (12) or higher
idx_core<- intersect(idx_core,which(nchar(gene_mhc1_mut$Peptide)+gene_mhc1_mut$Pos>10)) # Remove pos 0 (1) for 8/9/10-mers, ...
gene_mhc1_mut<- gene_mhc1_mut[idx_core,]
# Get ranks
gene_mhc1_mut_rank<- gene_mhc1_mut[,grep("Rank",colnames(gene_mhc1_mut))]
colnames(gene_mhc1_mut_rank)<- setdiff(unlist(strsplit(readLines(gene_mhc1_filename_mut,n=1),"\t")),"")
# Calculate PBR & PHBR
PBR_tmp<- apply(gene_mhc1_mut_rank,MARGIN = 2, function(x) tapply(x, gene_mhc1_mut$ID,"min"))
PHBR_proto_mhc1_mut<- apply(PBR_tmp,1,harmonic_mean)
# Sometimes peptides in beginning of protein and no 38 used for analysis, some by coincidence used multiple times; filter on only 38
pep_excl<- names(table(gene_mhc1_mut$ID)[table(gene_mhc1_mut$ID)!=38]) # 516
PHBR_proto_mhc1_mut[pep_excl]<- NA
# Add
GPPM_rand$PHBR_mut<- PHBR_proto_mhc1_mut[GPPM_rand$mut]

# dPHBR
GPPM_rand$dPHBR<- GPPM_rand$PHBR_mut-GPPM_rand$PHBR_wt

# Save
saveRDS(GPPM_rand,file="data/GPPM_rand.rds")
