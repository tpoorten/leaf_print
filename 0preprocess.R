library(plotly)
library(here)
source(paste0(here(),"/src/functions.R"))

## Read in data
dat = leaf_read(snp_data_file   = "input_data/SNP_data.txt", 
               sample_qc_file   = "input_data/Sample_QC.txt",
               # sample_meta_file = "input_data/sample_master_list.csv",
               sample_meta_file = "input_data/sampleInfo_all_cels_noQC.csv",
               ps_qc_file       = "input_data/Ps.performance.txt",
               snp_meta_file    = "input_data/snp_metadata_v2.txt"
              )

head(dat$snp_data[,1:3])
head(dat$ps_qc)
head(dat$sample_data)
View(dat$sample_data)
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


## Filter genotype calls with Confidence Scores
dat3 = leaf_confidence_filtering(leaf_data = dat2, 
                                 confidence_file = "input_data/confidences_filterRecommended.txt", 
                                 cutoff.confidence = 0.01)
View(dat3$sample_data)
dim(dat3$snp_data)
# //


## Create a gds file
library(SNPRelate)
# table(as.matrix(snp))
# length(unique(colnames(snp))) == ncol(snp)
rm(dat1, dat2)
table((dat3$ps_qc$istraw90_Chromsome))
table(is.na(dat3$ps_qc$istraw90_Chromsome))
table((dat3$ps_qc$bwaFv4.0_RNAME))
table(is.na(dat3$ps_qc$bwaFv4.0_RNAME))
#
dir.create("RDataset", showWarnings = F)
snpgdsCreateGeno("RDataset/snp_confidence.gds", genmat = as.matrix(dat3$snp_data),
                 sample.id = colnames(dat3$snp_data), snp.id = rownames(dat3$snp_data),
                 snp.chromosome = dat3$ps_qc$istraw90_Chromsome,
                 snp.position = dat3$ps_qc$istraw90_Position,
                 snp.allele = paste(dat3$ps_qc[,c("Allele.Ref")], dat3$ps_qc[,c("Allele.Alt")], sep = "/"), snpfirstdim=TRUE)
# make plink
# genofile <- snpgdsOpen("RDataset/snp.gds")
# geno.sample.ids = read.gdsn(index.gdsn(genofile, "sample.id"))
# snpgdsGDS2PED(genofile, "RDataset/snp", snp.id = dat3$ps_qc$probeset_id)
# snpgdsClose(genofile)
# //


## Analyze samples that match in dataset
# KING method (can incorporate family info) - MANUAL SAY USE ALL SNPS
genofile <- snpgdsOpen("RDataset/snp_confidence.gds", readonly = FALSE)

ibd.robust <- snpgdsIBDKING(genofile, num.thread=2, type=c("KING-robust"), autosome.only = F)
ibd.robust.coeff <- snpgdsIBDSelection(ibd.robust)
snpgdsClose(genofile)

ggplot(ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.02),], aes(IBS0, kinship)) + 
  geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed")
#
#


## Run analysis on sample matches: pedigree confirmation, ambiguous ID detection.
dat4 = leaf_analyze_matches(leaf_data = dat3, ibd.robust.coeff = ibd.robust.coeff, kinship.cutoff = 0.44, IBS0.cutoff = 0.002)
View(dat4$sample_data)
# //

## Create a gds file for ambiguous IDs
snpgdsCreateGeno("RDataset/snp_confidence_amb.gds", genmat = as.matrix(dat4$snp_data),
                 sample.id = colnames(dat4$snp_data), snp.id = rownames(dat4$snp_data),
                 snp.chromosome = dat4$ps_qc$istraw90_Chromsome,
                 snp.position = dat4$ps_qc$istraw90_Position,
                 snp.allele = paste(dat4$ps_qc[,c("Allele.Ref")], dat4$ps_qc[,c("Allele.Alt")], sep = "/"), snpfirstdim=TRUE)
# //

## Run trio pedigree analysis

# filter IDs - only include include IDs that appear to be closely related to another ID
genofile <- snpgdsOpen("RDataset/snp_confidence_amb.gds", readonly = FALSE)
geno.sample.ids = dat4$sample_data$Sample.Filename[which(dat4$sample_data$UseDownstream)]
  
ibd.robust <- snpgdsIBDKING(genofile, num.thread=2, type=c("KING-robust"), autosome.only = F, sample.id = geno.sample.ids)
ibd.robust.coeff <- snpgdsIBDSelection(ibd.robust)
snpgdsClose(genofile)

ggplot(ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.02),], aes(IBS0, kinship)) + 
  geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed")

ibd.robust.coeff.filter = ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.01),]
ids.filter = unique(c(ibd.robust.coeff.filter$ID1, ibd.robust.coeff.filter$ID2))
rm(ibd.robust.coeff)

# Get trio transgression ratios
trio.df = leaf_trio_pedigree(leaf_data = dat4, ibd.coeff = ibd.robust.coeff.filter, sample.ids = ids.filter)
trio.df = trio.df[which(trio.df$ttr < 0.03),]

# skip cases were Offspring and a Parent are replicates
trio.df = trio.df[ ifelse(gsub("_rep[A-Z]","", trio.df$ID.name) == gsub("_rep[A-Z]","", trio.df$Parent1.name) | 
                            gsub("_rep[A-Z]","", trio.df$ID.name) == gsub("_rep[A-Z]","", trio.df$Parent2.name), FALSE, TRUE)
                   ,]

# check Year info
trio.df$ID.year = dat4$sample_data$Year[match(trio.df$ID, dat4$sample_data$Sample.Filename)]
trio.df$check.year = ifelse({
  dat4$sample_data$Year[match(trio.df$ID, dat4$sample_data$Sample.Filename)] > 
    dat4$sample_data$Year[match(trio.df$Parent1, dat4$sample_data$Sample.Filename)] &
    dat4$sample_data$Year[match(trio.df$ID, dat4$sample_data$Sample.Filename)] > 
    dat4$sample_data$Year[match(trio.df$Parent2, dat4$sample_data$Sample.Filename)]
},"OK, offspring is\n more recent than parent", "NOT OK, parent(s) \nmore recent than offspring")
#
ggplot(trio.df,
       aes(ttr, total.tests, text=trio, color=check.year)) + geom_point()
ggplotly()

#
ggplot(trio.df,
       aes(ttr, total.tests, text=trio, color=ID.year)) + geom_point(size=1) + scale_color_gradient(low = "red", high="green")
ggplotly()

# //


## Run duo pedigree analysis

# filter IDs - only include include IDs that appear to be closely related to another ID
genofile <- snpgdsOpen("RDataset/snp_confidence_amb.gds", readonly = FALSE)
geno.sample.ids = dat4$sample_data$Sample.Filename[which(dat4$sample_data$UseDownstream)]

ibd.robust <- snpgdsIBDKING(genofile, num.thread=2, type=c("KING-robust"), autosome.only = F, sample.id = geno.sample.ids)
ibd.robust.coeff <- snpgdsIBDSelection(ibd.robust)
snpgdsClose(genofile)

ggplot(ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.02),], aes(IBS0, kinship)) + 
  geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed")

ibd.robust.coeff.filter = ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.01),]
ids.filter = unique(c(ibd.robust.coeff.filter$ID1, ibd.robust.coeff.filter$ID2))
rm(ibd.robust, ibd.robust.coeff)

# Process IBD-KING results
ibd.robust.coeff.filter.res = leaf_duo_ibd_king_process(leaf_data = dat4, ibd.coeff = ibd.robust.coeff.filter)

ggplot(ibd.robust.coeff.filter.res, aes(IBS0, kinship, color=relation, text=pair)) + geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed")
ggplotly()

ggplot(ibd.robust.coeff.filter.res, aes(IBS0, kinship, color=KingCall, text=pair)) + geom_point(alpha=0.7,size=.5) +
  geom_vline(xintercept = 0.001, linetype="dashed") + geom_vline(xintercept = 0.002, linetype="dashed")
ggplotly()
# //

# Get parent calls
ibd.robust.calls = leaf_duo_ibd_king_pedigree(leaf_data = dat4, ibd.coeff = ibd.robust.coeff.filter.res)

duo.res = ibd.robust.calls$duo.results.pairs
duo.ped = ibd.robust.calls$duo.results.pedigree

# //