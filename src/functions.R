read_leaf = function(snp_data_file=NULL, sample_qc_file=NULL, ps_qc_file = NULL, sample_meta_file=NULL){
  if(file.exists(snp_data_file)){
    snp_data = read.table(snp_data_file, sep="\t", header=T, stringsAsFactors = F, check.names = F)
    rownames(snp_data) = snp_data$probeset_id
    snp_data = snp_data[,which(colnames(snp_data) != "probeset_id")]
  } else {stop("Error: invalid filename for 'snp_data_file'")}
  
  if(file.exists(sample_qc_file)){
    sample_qc = read.table(sample_qc_file, sep="\t", header=T, stringsAsFactors = F)
    sample_qc = sample_qc[,1:7]
  } else {stop("Error: invalid filename for 'sample_qc_file'")}
  
  if(file.exists(ps_qc_file)){
    ps_qc = read.table(ps_qc_file, sep="\t", header=T, stringsAsFactors = F)
    rownames(ps_qc) = ps_qc$probeset_id
    if(identical(sort(rownames(snp_data)),sort(ps_qc$probeset_id))){
      ps_qc = ps_qc[match(rownames(snp_data),ps_qc$probeset_id),]
    } else {stop("Error: probeset IDs do not match between 'snp_data_file' and 'ps_qc_file'")}
  } else {stop("Error: invalid filename for 'ps_qc_file'")}
    
  sample_metadata = NULL
  if(!is.null(sample_meta_file)){
    if(file.exists(sample_meta_file)){
      sample_metadata = read.csv(sample_meta_file, stringsAsFactors = F)
      sample_metadata = sample_metadata[match(colnames(snp_data),sample_metadata[,1]),]
    } else {stop("Error: invalid filename for sample_meta_file. This file can be omitted.")}
  }
  
      
  leaf = list(snp_data = snp_data, sample_qc = sample_qc, ps_qc = ps_qc, sample_metadata = sample_metadata)
  return(leaf)
}


recalculate_sample_metrics = function(dataList = NULL){
  sample_qc = dataList$sample_qc
  snp_data = dataList$snp_data
  
  # Re-calculated call rate
  call_rate_recalculated = round(100 * apply(dataList$snp_data, 2, function(x) length(which(x >= 0))) / nrow(dataList$snp_data), digits = 3)
  dataList$sample_qc$call_rate_recalculated = call_rate_recalculated[match(dataList$sample_qc$Sample.Filename, names(call_rate_recalculated))]
  
  # Re-calculated heterozygous rate
  het_rate_recalculated = round(100 * apply(dataList$snp_data, 2, function(x) length(which(x == 1))) / apply(dataList$snp_data, 2, function(x) length(which(x >= 0))), digits = 3)
  dataList$sample_qc$het_rate_recalculated = het_rate_recalculated[match(dataList$sample_qc$Sample.Filename, names(het_rate_recalculated))]
  
  leaf = list(snp_data        = dataList$snp_data, 
              sample_qc       = sample_qc, 
              ps_qc           = dataList$ps_qc, 
              sample_metadata = dataList$sample_metadata)
  return(leaf)
}
  