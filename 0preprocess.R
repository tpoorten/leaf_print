rm(list=ls())
library(here)
source(paste0(here(),"/src/functions.R"))

## Read in data
dat = leaf_read(snp_data_file   = "input_data/SNP_data.txt", 
               sample_qc_file   = "input_data/Sample_QC.txt",
               ps_qc_file       = "input_data/Ps.performance.txt",
               sample_meta_file = "input_data/sample_master_list.csv"
              )

head(dat$snp_data)
head(dat$ps_qc)
head(dat$sample_data)
# //
