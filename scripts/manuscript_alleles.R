# Calculate OR for each allele? --> ALSO higher for HLA-C,

# Get drive rmutation
driver_muts<- colnames(readRDS("data/mut_matrix.rds"))# 

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

library(gplots)
muts<- readRDS("data/driver_muts13.rds")

pdf("temp/hm.pdf")
par(mfrow=c(1,3))
for(HLA in c("HLA-A","HLA-B","HLA-C")){
  data<- data.matrix(PBR)
  data<- data[,grep(HLA,colnames(data))]
  data<- na.omit(data)
  data<- data/apply(data,2,"median")
  data<- log2(data)
  data<- t(apply(data,1,"sort"))
  data<- data[muts,]
  data<- data[order(rowSums(data>0),decreasing = T),]
  data[data>3]<- 3
  data[data< -3]<- -3
  
  # hm<- heatmap.2(data, col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="column",density.info='none', distfun=function(x) dist(x,method = "manhattan"))
  # hm<- heatmap.2(data, col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none')
  hm<- heatmap.2(data, col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none',Rowv = F,Colv = F)
  # hm<- heatmap.2(data, col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="none",density.info='none', distfun=function(x) dist(x,method = "manhattan"))
  cat(rowMeans(data>0),"\n")
}
dev.off()

# Cor plot
data<- na.omit(PBR)
cor_matrix<- matrix(NA,ncol(data),ncol(data),dimnames = list(colnames(PBR),colnames(PBR)))
for(i in 1:nrow(cor_matrix)){
  for(j in 1:ncol(cor_matrix)){
    cor_tmp<- cor.test(data[,i],data[,j])
    cor_matrix[i,j]<- cor_tmp$estimate
  }
}
heatmap.2(cor_matrix,col=bluered(75),symkey=TRUE, key=TRUE, keysize=1, trace="none",scale="column",density.info='none', cexRow=0.5, cexCol=0.8,key.title = "Pearson's R",distfun=function(x) dist(x,method = "manhattan"))


# Rare alleles? --> Plug in 1kg
load("../immunoediting_2019/data/TCGA_HLA_types.RData")
freq_A<- sort(prop.table(table(c(TCGA_HLA_types_netMHC[,1],TCGA_HLA_types_netMHC[,2]))),decreasing = F)[1:2]
freq_B<- sort(prop.table(table(c(TCGA_HLA_types_netMHC[,3],TCGA_HLA_types_netMHC[,4]))),decreasing = F)[1:2]
freq_C<- sort(prop.table(table(c(TCGA_HLA_types_netMHC[,5],TCGA_HLA_types_netMHC[,6]))),decreasing = F)[1:2]
HLA_proto_mhc1<- c(names(freq_A),names(freq_B),names(freq_C))
# Same effect --> Argues against evolution


# 
vp_mut_matrix<- readRDS(file="data/vp_mut_matrix.rds")
mut_freq<- 100*colSums(vp_mut_matrix)/nrow(vp_mut_matrix)
cor.test(1:nrow(data),mut_freq[rownames(data)],method="spearman")
plot(mut_freq[rownames(data)])
# 
# hm$rowInd
# 
# clust_classes<- cutree(as.hclust(hm$rowDendrogram),4)
# table(clust_classes)
# sort(names(clust_classes[clust_classes==1]))
# 
# idx_proto<- which(colnames(PBR)%in%c("HLA-A01:01","HLA-A02:01","HLA-B07:02","HLA-B08:01","HLA-C07:01","HLA-C07:02"))
# rowMeans(PBR[1:20,idx_proto])/rowMeans(PBR[1:20,-idx_proto])
# 
# # Calculate OR all alleles?
# vp_mut_matrix<- readRDS("data/vp_mut_matrix.rds")
# vp_HLA1_mut_matrix<- vp_mut_matrix
# vp_HLA1_mut_matrix[,]<- NA
# common_muts<- intersect(rownames(PBR),colnames(vp_HLA1_mut_matrix))
# for(i in 1:nrow(vp_HLA1_mut_matrix)){
#   vp_HLA1_mut_matrix[i,common_muts]<- PBR[common_muts,1]
# }
# 
# vp_HLA_matrix_sel<- vp_HLA1_mut_matrix
# vp_mut_matrix_sel<- vp_mut_matrix[,colnames(vp_HLA_matrix_sel)]
# 
# source("scripts/functions/do_logreg.R")
# logreg_mut<- do_logreg(vp_mut_matrix_sel,vp_HLA_matrix_sel,"wp") # 1.09, takes +/- 0,5h
# tapply(as.numeric(vp_HLA_matrix_sel), as.numeric(vp_mut_matrix_sel), "median",na.rm=T) # 3,9 vs 5,0
# 
# # Turn around: PBR for 13 versus rest for all alleles, correlation frequency?
# driver_muts13 <- readRDS("~/projects/mhc_driver/data/driver_muts13.rds")
# idx_driver<- which(rownames(PBR)%in%driver_muts13)
# barplot(sort(colMeans(PBR[idx_driver,],na.rm=T)))
# barplot(sort(colMeans(PBR[-idx_driver,],na.rm=T)))
# 
# barplot(sort(log10(colMeans(PBR[idx_driver,],na.rm=T)/colMeans(PBR[-idx_driver,],na.rm=T)))) # Almost never lower :-)
# # Correlation frequency?


# # Load driver mutation data (needed to get ENSP id)
# load("temp/TCGA_maf_driver.RData")
# 
# # Harmonic mean function
# source("scripts/functions/harmonic_mean.R")
# 
# # Load TCGA Polysolver HLA genotypes & save all alleles to txt file to use with netMHCPan (max 1024 chars per line!)
# load("~/projects/immunoediting_2019/data/TCGA_HLA_types.RData")
# HLA_alleles_TCGA<- sort(setdiff(na.omit(unique(c(TCGA_HLA_types_netMHC[,1],TCGA_HLA_types_netMHC[,2],TCGA_HLA_types_netMHC[,3],TCGA_HLA_types_netMHC[,4],TCGA_HLA_types_netMHC[,5],TCGA_HLA_types_netMHC[,6]))),""))
# HLA_alleles_TCGA<- setdiff(HLA_alleles_TCGA,"HLA-B13:07") # B13:07 not in netMHCPan!
# cat(HLA_alleles_TCGA[1:65],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt")
# cat("\n",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
# cat(HLA_alleles_TCGA[66:130],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
# cat("\n",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
# cat(HLA_alleles_TCGA[131:length(HLA_alleles_TCGA)],sep=",",file = "temp/netMHCPan40/TCGA_HLA_alleles.txt",append = T)
# 
# # Create empty HLA matrix
# HLA_matrix<- mut_matrix
# HLA_matrix[,]<- NA
# 
# # Get peptides for driver muts
# library(EnsDb.Hsapiens.v86) # To get ensembl protein sequence information
# edb <- filter(EnsDb.Hsapiens.v86, filter = AnnotationFilterList(TxBiotypeFilter("protein_coding"))) # Restrict edb to coding genes
# for(mut_tmp in colnames(HLA_matrix)[1:ncol(HLA_matrix)]){
#   
#   cat(grep(mut_tmp,colnames(HLA_matrix)),"\n")
#   mut_tmp_gene<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Hugo_Symbol"])
#   mut_tmp_ENSP<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"ENSP"])
#   mut_tmp_aa<- unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Amino_acids"])
#   mut_tmp_aa_ref<- gsub("\\/.","",mut_tmp_aa)
#   mut_tmp_aa_alt<- gsub(".\\/","",mut_tmp_aa)
#   mut_tmp_pos<- as.numeric(unique(TCGA_maf_driver[TCGA_maf_driver$mut_id==mut_tmp,"Protein_position"]))
#   
#   # Get protein sequence information
#   if(mut_tmp_gene=="RUNX1T1") mut_tmp_ENSP<- "ENSP00000428742" # Manual correction, mutations only matched this ENSP
#   if(mut_tmp_gene=="CRLF2") mut_tmp_ENSP<- "ENSP00000370978" # Manual correction, mutations only matched this ENSP
#   if(mut_tmp_gene=="TCF7L2") mut_tmp_ENSP<- "ENSP00000486891" # Manual correction, mutations only matched this ENSP
#   if(mut_tmp_gene=="MYC") mut_tmp_ENSP<- "ENSP00000478887" # Manual correction, mutations only matched this ENSP
#   if(mut_tmp_ENSP=="ENSP00000298047") next # not found
#   prot_sequence<- proteins(edb, filter = ProteinIdFilter(mut_tmp_ENSP), return.type = "AAStringSet")
#   # Use gene if ENSP not found
#   if(length(prot_sequence)==0){ 
#     proteins_for_gene<- proteins(edb, filter = GeneNameFilter(mut_tmp_gene), return.type = "AAStringSet")
#     prot_sequence <- proteins_for_gene[which(width(proteins_for_gene)==max(width(proteins_for_gene)))[1]] # Take max if multiple
#   }
#   
#   # Get 21-mers
#   sequence_strings <- as.character(prot_sequence)
#   pep11_tmp<- substr(sequence_strings,mut_tmp_pos-10,mut_tmp_pos+10)
#   # IGV check: ok for BRAF V600E
#   
#   # Add mutations
#   pos_in_pep_tmp<- 11
#   if(mut_tmp_pos<11) pos_in_pep_tmp<-  mut_tmp_pos
#   if(substr(pep11_tmp,pos_in_pep_tmp,pos_in_pep_tmp)!=mut_tmp_aa_ref) stop("Check mutation: ", mut_tmp)
#   # Adapt the following after errors: "APCS2307L"
#   substr(pep11_tmp,pos_in_pep_tmp,pos_in_pep_tmp)<- mut_tmp_aa_alt
#   
#   # Save peptide faste files
#   cat(">",mut_tmp,"\n",pep11_tmp,"\n",sep="",file = paste0("temp/netMHCPan40/input/",mut_tmp,".fasta"))
#   
#   # Run NetMHCPan: split in 3 because netMHCPan only accepts HLA allele file with fewer than 1024 chars
#   system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -1 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"1.xls > /dev/null"))
#   system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -2 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"2.xls > /dev/null"))
#   system(paste0("/home/labgroups/ccgg/tools/netMHCpan-4.0/netMHCpan -a `head -3 temp/netMHCPan40/TCGA_HLA_alleles.txt | tail -1` -f temp/netMHCPan40/input/",mut_tmp,".fasta  -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile temp/netMHCPan40/output/",mut_tmp,"3.xls > /dev/null"))
#   
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
#   HLA_matrix[common_samples,mut_tmp]<- PHBR_tmp[common_samples]
# }
# 
# # Save
# saveRDS(HLA_matrix,file="data/HLA_matrix.rds")

