allele_PBR_ls<- readRDS("results/data/TCGA_alleles_OR.rds")

# Data frame of results
allele_PBR_df<- as.data.frame(t(sapply(allele_PBR_ls,function(x) x)))

# add FDR
allele_PBR_df$q<- p.adjust(allele_PBR_df$p,"fdr")

# Get % OR increasing
mean(allele_PBR_df$q<0.05&allele_PBR_df$OR>1) # 84.5%
mean(allele_PBR_df$q<0.05&allele_PBR_df$OR<1) # 8.8%

# Get % OR by allele
allele_PBR_df$gene<- NA
allele_PBR_df$gene[grep("HLA-A",rownames(allele_PBR_df))]<- "HLA-A"
allele_PBR_df$gene[grep("HLA-B",rownames(allele_PBR_df))]<- "HLA-B"
allele_PBR_df$gene[grep("HLA-C",rownames(allele_PBR_df))]<- "HLA-C"

tapply(allele_PBR_df$q<0.05&allele_PBR_df$OR>1,allele_PBR_df$gene,"mean") 
# HLA-A     HLA-B     HLA-C 
# 0.6792453 0.8737864 1.0000000
tapply(allele_PBR_df$q<0.05&allele_PBR_df$OR>1,allele_PBR_df$gene,"sum") 
# HLA-A HLA-B HLA-C 
# 36    90    37 
tapply(allele_PBR_df$q<0.05&allele_PBR_df$OR>1,allele_PBR_df$gene,"length") 
# HLA-A HLA-B HLA-C 
# 53   103    37 

tapply(allele_PBR_df$q<0.05&allele_PBR_df$OR<1,allele_PBR_df$gene,"mean") 
# HLA-A      HLA-B      HLA-C 
# 0.18867925 0.06796117 0.00000000 
tapply(allele_PBR_df$q<0.05&allele_PBR_df$OR<1,allele_PBR_df$gene,"sum") 
# HLA-A HLA-B HLA-C 
# 10     7     0 

# Plot
allele_PBR_df<- allele_PBR_df[order(allele_PBR_df$OR),]
allele_PBR_df$gene<- factor(allele_PBR_df$gene)

svglite::svglite("results/figs/manuscript_alleles_OR.svg")
plot(allele_PBR_df$OR,frame.plot=F,axes=F,ylab="OR",xlab=NA,pch=16,col=as.numeric(allele_PBR_df$gene)+1,cex=0.8,ylim=c(0.9,1.25))
axis(2)
abline(h=1, lty=2)
for(i in 1:nrow(allele_PBR_df)) lines(c(i,i),c(allele_PBR_df$`2.5 %`[i],allele_PBR_df$`97.5 %`[i]),xpd=NA)
legend("topleft",levels(allele_PBR_df$gene),bty="n",fill=2:4)
dev.off()

# # correlation with population frequency?
# pop_HLA<- read.table("downloads/1kg/1000_genomes_hla.tsv",sep="\t",header = T,row.names = 3,stringsAsFactors = F)
# idx_amb<- unique(as.numeric(unlist(apply(pop_HLA[,3:8],2,function(x) grep("\\/",x))))) # identify ambiguous samples (e.g. 02/03 format): 63 in total
# pop_HLA<- pop_HLA[-idx_amb,] 
# pops<- pop_HLA$Region
# pop_HLA<- pop_HLA[,3:8] # Focus on MHC-I only (MHC-II not complete)
# freq_A<- sort(prop.table(table(c(pop_HLA$HLA.A.1,pop_HLA$HLA.A.2))),decreasing = T)
# names(freq_A)<- paste0("HLA-A",names(freq_A))
# freq_B<- sort(prop.table(table(c(pop_HLA$HLA.B.1,pop_HLA$HLA.B.2))),decreasing = T)
# names(freq_B)<- paste0("HLA-B",names(freq_B))
# freq_C<- sort(prop.table(table(c(pop_HLA$HLA.C.1,pop_HLA$HLA.C.2))),decreasing = T)
# names(freq_C)<- paste0("HLA-C",names(freq_C))
# freq<- c(freq_A,freq_B,freq_C)
# allele_PBR_df$freq<- freq[rownames(allele_PBR_df)]
# plot(allele_PBR_df$freq, allele_PBR_df$OR)
# cor.test(allele_PBR_df$freq, allele_PBR_df$OR)

######################################
# Driver muts alleles mut vs not mut
######################################

load("results/data/TCGA_alleles_PBR_t.RData")
library(ggpubr)

p_mut_PBR<- ggplot(PBR_t, aes(x=mut, y=PBR, fill=isMut)) +
  geom_boxplot(outlier.shape=NA) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=PBR_med_all) 

p_mut_PHBR<- ggplot(PHBR_t, aes(x=mut, y=PHBR, fill=isMut)) +
  geom_boxplot(outlier.shape=NA) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=PHBR_med_all)  + 
  stat_compare_means(method = "wilcox.test", paired = FALSE, aes(group = isMut),label = "p.format")

p_HLA<- ggplot(PBR_t, aes(x=mut, y=PBR, fill=HLA)) +
  geom_boxplot(outlier.shape=NA) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept=c(PBR_med_A,PBR_med_B,PBR_med_C))

ggsave("results/figs/TCGA_PHBR_driverMuts.svg",p_mut_PHBR)

