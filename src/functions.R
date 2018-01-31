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
