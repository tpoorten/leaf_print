## leaf_print

In early development phase. 

Developers: Tom Poorten, Randi Famula, Mitchell Feldmann

Functions to build:

* read in data
* filter SNPs - keep only polymorphic well-behaving (biallelic) SNPs
* make QAQC plots - sample-wise missing data distribution, ...?
* filter samples - by quality or vector of arbitrary sample IDs
* calc sample heterozygosity
* convert to GDS format for [SNPRelate](http://bioconductor.org/packages/release/bioc/html/SNPRelate.html) R package
* assess replicate samples and keep one sample per replicate set
  * for expected reps - keep sample with best call rate
  * for unexpected reps - keep sample that is pedigree-confirmed if possible
* pop structure viz - PCA
  * flag unexpected PC values based on reported taxonomy
* pedigree analysis (IBD-based inference of PO relationships)

