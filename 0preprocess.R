library(here)
source(paste0(here(),"/src/functions.R"))

## Read in data
dat = leaf_read(snp_data_file   = "input_data/SNP_data.txt", 
               sample_qc_file   = "input_data/Sample_QC.txt",
               ps_qc_file       = "input_data/Ps.performance.txt",
               sample_meta_file = "input_data/sample_master_list.csv"
              )

head(dat$snp_data)
nrow(dat$ps_qc)
head(dat$sample_data)
# //
dat1 = leaf_filter(leaf_data = dat)

head(dat1$snp_data)
head(dat1$ps_qc)
head(dat1$sample_data)

