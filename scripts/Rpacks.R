# Rpacks
#########

# Install additional R packages (from console, RStudio crashes in conda when installing packages)
install.packages("plotrix") # 2D matplot
install.packages("lme4") # 2D matplot
install.packages("gplots") # heatmap
install.packages("vioplot") # heatmap
install.packages("doParallel")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(ggpubr) # Function compare_stat_means
