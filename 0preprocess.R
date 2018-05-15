library(here)
source(paste0(here(),"/src/functions.R"))

## Read in data
dat = leaf_read(snp_data_file   = "input_data/SNP_data.txt", 
               sample_qc_file   = "input_data/Sample_QC.txt",
               # sample_meta_file = "input_data/sample_master_list.csv",
               sample_meta_file = "input_data/sampleInfo_all_cels.csv",
               ps_qc_file       = "input_data/Ps.performance.txt",
               snp_meta_file    = "input_data/snp_metadata_v2.txt"
              )

head(dat$snp_data[,1:3])
head(dat$ps_qc)
head(dat$sample_data)

# //

## Filter SNPs
dat1 = leaf_filter(leaf_data = dat, filter_col1 = "ConversionType", condition1 = "PolyHighResolution", filter_col2 = "BestProbeset", condition2 = 1)

dim(dat$snp_data)
dim(dat1$snp_data)
dim(dat$ps_qc)
dim(dat1$ps_qc)
head(dat1$sample_data)
# //

## Recalculate sample metrics
dat2 = leaf_recalculate_sample_metrics(leaf_data = dat1)

head(dat1$sample_data)
head(dat2$sample_data)
# //

## Make ID_reps columns with leaf_mark_replicate_IDs()
dat2 = leaf_mark_replicate_IDs(leaf_data = dat2)

head(dat2$sample_data)
View(dat2$sample_data)
# //


## Create a gds file
library(SNPRelate)
# table(as.matrix(snp))
# length(unique(colnames(snp))) == ncol(snp)
table((dat1$ps_qc$istraw90_Chromsome))
table(is.na(dat1$ps_qc$istraw90_Chromsome))
table((dat1$ps_qc$bwaFv4.0_RNAME))
table(is.na(dat1$ps_qc$bwaFv4.0_RNAME))
#
dir.create("RDataset", showWarnings = F)
snpgdsCreateGeno("RDataset/snp.gds", genmat = as.matrix(dat1$snp_data),
                 sample.id = colnames(dat1$snp_data), snp.id = rownames(dat1$snp_data),
                 snp.chromosome = dat1$ps_qc$istraw90_Chromsome,
                 snp.position = dat1$ps_qc$istraw90_Position,
                 snp.allele = paste(dat1$ps_qc[,c("Allele.Ref")], dat1$ps_qc[,c("Allele.Alt")], sep = "/"), snpfirstdim=TRUE)
# make plink
# genofile <- snpgdsOpen("RDataset/snp.gds")
# geno.sample.ids = read.gdsn(index.gdsn(genofile, "sample.id"))
# snpgdsGDS2PED(genofile, "RDataset/snp", snp.id = dat1$ps_qc$probeset_id)
# snpgdsClose(genofile)
