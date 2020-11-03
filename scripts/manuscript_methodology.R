pep_tmp<- rep(0,21)
# 8-mers
for(j in 4:11){
  pep_tmp[j:(j+7)]<- pep_tmp[j:(j+7)]+1
}
# 9-mers
for(j in 3:11){
  pep_tmp[j:(j+8)]<- pep_tmp[j:(j+8)]+1
}
# 10-mers
for(j in 2:11){
  pep_tmp[j:(j+9)]<- pep_tmp[j:(j+9)]+1
}
# 11-mers
for(j in 1:11){
  pep_tmp[j:(j+10)]<- pep_tmp[j:(j+10)]+1
}

# plot
svglite::svglite("results/figs/manuscript_methodology.svg")
# par(mfrow=c(2,1))
plot(pep_tmp,type="b",axes=F,xlab=NA,ylab="# peptides containing aa",ylim=c(0,80))
axis(1,at=1:21,1:21)
axis(2)

# Add peptides
n<- 80
points(1:21,y=rep(n,21))

# 8-mers
for(j in 4:11){
  n<- n-1
  points(j:(j+7),y=rep(n,8))
}
# 9-mers
for(j in 3:11){
  n<- n-1
  points(j:(j+8),y=rep(n,9))
}
# 10-mers
for(j in 2:11){
  n<- n-1
  points(j:(j+9),y=rep(n,10))
}
# 11-mers
for(j in 1:11){
  n<- n-1
  points(j:(j+10),y=rep(n,11))
}

dev.off()

# BRAF example
###############

mut_tmp<- "BRAF_V600E"

# MHC-I
for(i in c(1:3)){
  gene_mhc1_filename<- paste0("temp/netMHCPan40/output/",mut_tmp,i,'.xls')
  if(!file.exists(gene_mhc1_filename)) next
  gene_mhc1<- read.table(gene_mhc1_filename,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
  # Only select peptides with actual mutation present
  idx_core<- which(gene_mhc1$Pos%in%0:10) # Cannever start at position 11 (12) or higher
  idx_core<- intersect(idx_core,which(nchar(gene_mhc1$Peptide)+gene_mhc1$Pos>10)) # Remove pos 0 (1) for 8/9/10-mers, ...
  gene_mhc1<- gene_mhc1[idx_core,]
  # Get ranks
  gene_mhc1_rank_tmp<- gene_mhc1[,grep("Rank",colnames(gene_mhc1))]
  colnames(gene_mhc1_rank_tmp)<- setdiff(unlist(strsplit(readLines(gene_mhc1_filename,n=1),"\t")),"")
  # Calculate PBR
  if(i==1) gene_mhc1_rank<- gene_mhc1_rank_tmp
  else gene_mhc1_rank<- cbind(gene_mhc1_rank,gene_mhc1_rank_tmp)
}

# PR
HLA_proto_mhc1<- as.character(read.table("../mhc2/data/mhc1_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))
PR<- gene_mhc1_rank[,HLA_proto_mhc1]
PR_n<- nchar(gene_mhc1$Peptide)
PR_pos<- as.numeric(gene_mhc1$Pos)
PR<- PR[order(PR_n,PR_pos),]

# PBR
PBR<- apply(PR,MARGIN = 2, "min")
PBR
# HLA-A01:01 HLA-A02:01 HLA-B07:02 HLA-B08:01 HLA-C07:01 HLA-C07:02 
# 7.0888    25.1783    20.3467    15.1122     5.8338     6.2618 

# PHBR
source("scripts/functions/harmonic_mean.R")
PHBR<- harmonic_mean(PBR) # 9.566061 for prototypical (1.91 for wt)

# Save in excel
PR<- round(PR,1)
WriteXLS::WriteXLS("PR","results/tables/PR_BRAF.xlsx")

# MHC-II
gene_mhc2_filename<- paste0("temp/netmhciipan32//output/driver_mut.xls")
gene_mhc2<- read.table(gene_mhc2_filename,header = T, stringsAsFactors = F,sep ="\t",quote = "",skip=1)
gene_mhc2<- gene_mhc2[gene_mhc2$ID==mut_tmp,]
# Get ranks
gene_mhc2_rank_tmp<- gene_mhc2[,grep("Rank",colnames(gene_mhc2))]
colnames(gene_mhc2_rank_tmp)<- setdiff(unlist(strsplit(readLines(gene_mhc2_filename,n=1),"\t")),"")
# Calculate PBR
gene_mhc2_rank<- gene_mhc2_rank_tmp

# PR
HLA_proto_mhc2<- as.character(read.table("../mhc2/data/mhc2_HLA_alleles_proto.txt",sep=",",stringsAsFactors = F))
PR<- gene_mhc2_rank[,HLA_proto_mhc2]
PR_n<- nchar(gene_mhc2$Peptide)
PR_pos<- as.numeric(gene_mhc2$Pos)
PR<- PR[order(PR_n,PR_pos),]

# PBR
PBR<- apply(PR,MARGIN = 2, "min")
PBR
# HLA-DPA10103-DPB10201 HLA-DPA10103-DPB10401 HLA-DPA10201-DPB10201 HLA-DPA10201-DPB10401 
# 35                    39                    35                    39 
# HLA-DQA10102-DQB10202 HLA-DQA10102-DQB10301 HLA-DQA10501-DQB10202 HLA-DQA10501-DQB10301 
# 48                    36                    35                    28 
# DRB1_0701             DRB1_1501 
# 43                    50 

# PHBR
source("scripts/functions/harmonic_mean.R")
PHBR<- harmonic_mean(PBR) # 37.8

# Save in excel
PR<- round(PR,1)
WriteXLS::WriteXLS("PR","results/tables/PR_BRAF_MHCII.xlsx")