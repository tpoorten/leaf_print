#Functions for the Leaf_print project
#Tom Poorten, Mitchell Feldmann, Randi Famula

################################################
#The purpose of this function is to read and compile all sample, SNP, and metadata.
leaf_read = function(snp_data_file=NULL, ps_qc_file = NULL, sample_qc_file=NULL, sample_meta_file=NULL, snp_meta_file = NULL){
  print("Reading in SNP data ...")
  # SNP DATA
  if(file.exists(snp_data_file)){
    snp_data = read.table(snp_data_file, sep="\t", header=T, stringsAsFactors = F, check.names = F)
    colnames(snp_data) = gsub("_call_code","",colnames(snp_data))
    rownames(snp_data) = snp_data$probeset_id
    snp_data = snp_data[,which(colnames(snp_data) != "probeset_id")]
  } else {stop("Error: invalid filename for 'snp_data_file'")}
  
  print("Reading in metadata files ...")
  # SAMPLE_QC - from Axiom Analysis Suite output
  if(file.exists(sample_qc_file)){
    sample_qc = read.table(sample_qc_file, sep="\t", header=T, stringsAsFactors = F)
    sample_qc = sample_qc[,1:7]
    rownames(sample_qc) = sample_qc$Sample.Filename
    # keep on CELs in 'snp_data'
    sample_qc = sample_qc[match(colnames(snp_data), rownames(sample_qc)),]
  } else {stop("Error: invalid filename for 'sample_qc_file'")}
  
  # SAMPLE METADATA - any extra sample metadata to add to sample_qc
  sample_metadata = NULL
  # sample_all_data = NULL
  if(!is.null(sample_meta_file)){
    if(file.exists(sample_meta_file)){
      sample_metadata = read.csv(sample_meta_file, stringsAsFactors = F)
      sample_metadata = sample_metadata[match(colnames(snp_data),sample_metadata[,1]),]
      rownames(sample_metadata) = sample_metadata$Best.Array
      sample_qc = merge(sample_qc, sample_metadata, by = 1) # merge by column 1
      rownames(sample_qc) = sample_qc$Sample.Filename
    } else {stop("Error: invalid filename for sample_meta_file. This file can be omitted.")}
  }
  
  # PROBESET QC - from Axiom Analysis Suite output
  if(file.exists(ps_qc_file)){
    ps_qc = read.table(ps_qc_file, sep="\t", header=T, stringsAsFactors = F)
    rownames(ps_qc) = ps_qc$probeset_id
    if(identical(sort(rownames(snp_data)),sort(rownames(ps_qc)))){
      ps_qc = ps_qc[match(rownames(snp_data),rownames(ps_qc)),]
    } else {stop("Error: probeset IDs do not match between 'snp_data_file' and 'ps_qc_file'")}
  } else {stop("Error: invalid filename for 'ps_qc_file'")}
  
  # PROBESET METADATA - any extra snp metadata to add to sample_qc
  if(!is.null(snp_meta_file)){
    if(file.exists(snp_meta_file)){
      snp_metadata = read.table(snp_meta_file, header=T, sep="\t", stringsAsFactors = F)
      # rownames(snp_metadata) = snp_metadata$probeset_id
      ps_qc = merge(ps_qc, snp_metadata, by = 1) # merge by column 1
      rownames(ps_qc) = ps_qc$probeset_id
    } else {stop("Error: invalid filename for sample_meta_file. This file can be omitted.")}
  }
  
  print(paste0("snp_data   (data.frame)  : ",nrow(snp_data)," SNPs (rows)  x  ", ncol(snp_data), " CELs (columns)"))
  print(paste0("ps_qc      (data.frame)  : ",nrow(ps_qc)," SNPs (rows)  x  ", ncol(ps_qc), " info columns"))
  print(paste0("sample_qc  (data.frame)  : ",nrow(sample_qc)," samples (rows)  x  ", ncol(sample_qc), " info columns"))

  if(identical(rownames(ps_qc), rownames(snp_data))){
    print("    Rows in 'ps_qc' matches up with rows in 'snp_data'")
  } else {
    print("    Rows in 'ps_qc' DOES NOT match up with rows in 'snp_data'")
  }
  
  if(identical(rownames(sample_qc), colnames(snp_data))){
    print("    Rows in 'sample_qc' matches up with columns in 'snp_data'")
  } else {
    print("    Rows in 'sample_qc' DOES NOT match up with columns in 'snp_data'")
  }
  print("DONE")
  
  leaf = list(snp_data        = snp_data, 
              ps_qc           = ps_qc, 
              sample_data     = sample_qc)
  
  return(leaf)
}

########################################
# the point of this function is to filter snp_data and ps_qc by values in ConversionType and BestProbeset in ps_qc
# this fxn will filter on the columns specified in filter_col1 and filter_col2 by condition1 and condition2
leaf_filter = function(leaf_data = NULL, filter_col1 = "ConversionType", condition1 = "PolyHighResolution", filter_col2 = "BestProbeset", condition2 = 1){
  
  if(!is.null(leaf_data)){
    
    if(!is.null(filter_col1)){
    leaf_data$ps_qc = leaf_data$ps_qc[which(leaf_data$ps_qc[,filter_col1] == condition1),]
    leaf_data$snp_data = leaf_data$snp_data[match(rownames(leaf_data$ps_qc), rownames(leaf_data$snp_data)),]
    } else {stop("Error: Require at least one column and one condition")}
    
    if(!is.null(filter_col2)){
      leaf_data$ps_qc = leaf_data$ps_qc[which(leaf_data$ps_qc[,filter_col2] == condition2),]
      leaf_data$snp_data = leaf_data$snp_data[match(rownames(leaf_data$ps_qc), rownames(leaf_data$snp_data)),]
    } else {}
    
  } else {stop("Error: Object Not Found")}
  
  if(identical(rownames(leaf_data$ps_qc), rownames(leaf_data$snp_data))){
    print("    Rows in 'ps_qc' matches up with rows in 'snp_data'")
  } else {
    print("    Rows in 'ps_qc' DOES NOT match up with rows in 'snp_data'")
  }
  
  if(identical(rownames(leaf_data$sample_data), colnames(leaf_data$snp_data))){
    print("    Rows in 'sample_qc' matches up with columns in 'snp_data'")
  } else {
    print("    Rows in 'sample_qc' DOES NOT match up with columns in 'snp_data'")
  }
  
  leaf = list(snp_data        = leaf_data$snp_data, 
              ps_qc           = leaf_data$ps_qc, 
              sample_data     = leaf_data$sample_data)
  return(leaf)
}

########################################
#Purpose is to recalculate sample call and heterozygous call rates

leaf_recalculate_sample_metrics = function(dataList = NULL){
  
  # Re-calculated call rate
  call_rate_recalculated = round(100 * apply(dataList$snp_data, 2, function(x) length(which(x >= 0))) / nrow(dataList$snp_data), digits = 3)
  dataList$sample_data$call_rate_recalculated = call_rate_recalculated[match(dataList$sample_data$Sample.Filename, names(call_rate_recalculated))]
  
  # Re-calculated heterozygous rate
  het_rate_recalculated = round(100 * apply(dataList$snp_data, 2, function(x) length(which(x == 1))) / apply(dataList$snp_data, 2, function(x) length(which(x >= 0))), digits = 3)
  dataList$sample_data$het_rate_recalculated = het_rate_recalculated[match(dataList$sample_data$Sample.Filename, names(het_rate_recalculated))]
  
  return(dataList)
}
  

