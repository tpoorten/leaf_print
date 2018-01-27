read_leaf = function(snp_data_file="", sample_qc_file="", ps_qc_file = "", sample_meta_file=""){
  snp_data = read.table(snp_data_file, sep="\t", header=T, stringsAsFactors = F, check.names = F)
  rownames(snp_data) = snp_data$probeset_id
  snp_data = snp_data[,which(colnames(snp_data) != "probeset_id")]
  
  sample_metadata = read.csv(sample_meta_file, stringsAsFactors = F)
  sample_metadata = sample_metadata[match(colnames(snp_data),sample_metadata[,1]),]
  
  sample_qc = read.table(sample_qc_file, sep="\t", header=T, stringsAsFactors = F)
  sample_qc = sample_qc[,1:7]
  
  ps_qc = read.table(ps_qc_file, sep="\t", header=T, stringsAsFactors = F)
  rownames(ps_qc) = ps_qc$probeset_id
  ps_qc = ps_qc[match(rownames(snp_data),ps_qc$probeset_id),]
  
  
  leaf = list(snp_data = snp_data, sample_qc = sample_qc, ps_qc = ps_qc, sample_metadata = sample_metadata)
  return(leaf)
}
