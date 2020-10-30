# Get all download scripts
# ln downloads/TCGA/maf/mc3/download_process.R scripts/downloads/download_mc3.R # Hard links don't work, different partiion

# Just copy
cp downloads/TCGA/maf/mc3/download_process.R scripts/downloads/download_mc3.R 
cp downloads/cosmic/CGC/v91/download_process.R scripts/downloads/download_cgc.R
cp downloads/TCGA/clin/download_process.R scripts/downloads/download_TCGA_id.R 
