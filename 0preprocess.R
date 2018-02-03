library(here)
source(paste0(here(),"/src/functions.R"))

## Read in data
dat = leaf_read(snp_data_file   = "input_data/SNP_data.txt", 
               sample_qc_file   = "input_data/Sample_QC.txt",
               ps_qc_file       = "input_data/Ps.performance.txt",
               sample_meta_file = "input_data/sample_master_list.csv"
              )

head(dat$snp_data)
table(dat$ps_qc$H.W.p.Value)
head(dat$sample_data)

# //

dat1 = leaf_filter(leaf_data = dat, filter_col1 = "ConversionType", condition1 = "PolyHighResolution", filter_col2 = "BestProbeset", condition2 = 1)

head(dat1$snp_data)
nrow(dat1$ps_qc)
head(dat1$sample_data)

## 
dat2 = recalculate_sample_metrics(dataList = dat)

