#Functions for the Leaf_print project
#Tom Poorten, Mitchell Feldmann, Randi Famula

# To Do
# confid - make function to make figure of confidence scores
# leaf_analyze_matches() - knownPO or poKnown?

################################################
# The purpose of this function is to read and compile all sample, SNP, and metadata.
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
# Purpose is to recalculate sample call and heterozygous call rates

leaf_recalculate_sample_metrics = function(leaf_data = NULL){
  print("    Adding new columns to leaf$sample_data: call_rate_recalculated, het_rate_recalculated")
  
  # Re-calculated call rate
  call_rate_recalculated = round(100 * apply(leaf_data$snp_data, 2, function(x) length(which(x >= 0))) / nrow(leaf_data$snp_data), digits = 3)
  leaf_data$sample_data$call_rate_recalculated = call_rate_recalculated[match(leaf_data$sample_data$Sample.Filename, names(call_rate_recalculated))]
  
  # Re-calculated heterozygous rate
  het_rate_recalculated = round(100 * apply(leaf_data$snp_data, 2, function(x) length(which(x == 1))) / apply(leaf_data$snp_data, 2, function(x) length(which(x >= 0))), digits = 3)
  leaf_data$sample_data$het_rate_recalculated = het_rate_recalculated[match(leaf_data$sample_data$Sample.Filename, names(het_rate_recalculated))]
  
  return(leaf_data)
}
  

########################################
# Automatically add on suffix for replicates
leaf_mark_replicate_IDs = function(leaf_data = NULL, id_column_name = "ID"){
  print("    Adding new column to leaf$sample_data: ID_reps")
  leaf_data$sample_data = leaf_data$sample_data[order(leaf_data$sample_data$ID),]
  leaf_data$sample_data$ID_reps = leaf_data$sample_data$ID
  leaf_data$sample_data$rep = FALSE
  leaf_data$sample_data$rep[which(leaf_data$sample_data$ID %in% names(table(leaf_data$sample_data$ID))[which(table(leaf_data$sample_data$ID) > 1)])] = TRUE
  leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$rep)] = unlist(tapply(leaf_data$sample_data$ID[which(leaf_data$sample_data$rep)], leaf_data$sample_data$ID[which(leaf_data$sample_data$rep)], function(x) paste(x,LETTERS[1:length(x)], sep = "_rep")))

  if(! identical(rownames(leaf_data$sample_data), colnames(leaf_data$snp_data))){
    print("    Re-sorting columns in leaf$snp_data based on sample order in leaf$sample_data")
    leaf_data$snp_data = leaf_data$snp_data[,match(rownames(leaf_data$sample_data), colnames(leaf_data$snp_data))]
  }

  return(leaf_data)
}
  
#########################################
# Use SNP confidence - turn low confidence to missing data
leaf_confidence_filtering = function(leaf_data = NULL, confidence_file = NULL, cutoff.confidence = 0.01){
  print("    Reading in Confidence Scores")
  snpCondfid = read.table(confidence_file, header=T, check.names = F, stringsAsFactors = F)
  
  rownames(snpCondfid) = snpCondfid$probeset_id
  snpCondfid = snpCondfid[,grep("confidence", colnames(snpCondfid))]
  colnames(snpCondfid) = sapply(colnames(snpCondfid), function(x) unlist(strsplit(x, "_confid"))[1])

  ## make snpConfid match snp matrix
  snpsKeep = intersect(rownames(leaf_data$snp_data), rownames(snpCondfid))
  snpCondfid = snpCondfid[which(rownames(snpCondfid) %in% snpsKeep),]
  snp.confidFilter = leaf_data$snp_data[which(rownames(leaf_data$snp_data) %in% snpsKeep),]
  snpCondfid = (snpCondfid[match(rownames(snp.confidFilter), rownames(snpCondfid)),
                           match(colnames(snp.confidFilter), colnames(snpCondfid))] )

  # QC
  snpCondfidColMean = colMeans(snpCondfid, na.rm = T)
  leaf_data$sample_data$MeanSNPconfid = snpCondfidColMean[match(leaf_data$sample_data$Sample.Filename, names(snpCondfidColMean))]

  hist(snpCondfid[snpCondfid > 0.005], xlim=c(0,.2), breaks=200)
  
  # length(snpCondfid[snpCondfid > 0.005])

  ## make genotype changes
  if(identical(rownames(snp.confidFilter), rownames(snpCondfid)) & identical(colnames(snp.confidFilter), colnames(snpCondfid))){
    print("Data Frames Match; Proceed with filtering")
    snp.confidFilter[snpCondfid >= cutoff.confidence] = -1
  } else{print("ERROR with Data Frame Matching!")}
  
  leaf_data$snp_data = snp.confidFilter
  leaf_data$ps_qc = leaf_data$ps_qc[match(rownames(leaf_data$snp_data), rownames(leaf_data$ps_qc)),]
  
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
  
  return(leaf_data)
}





#########################################
# Analyze matches
leaf_analyze_matches = function(leaf_data = NULL, ibd.robust.coeff = NULL, kinship.cutoff = 0.44, IBS0.cutoff = 0.002, matches.samples.keep = NULL, matches.samples.exclude = NULL){
  # leaf_data = dat1
  # kinship.cutoff = 0.44
  # IBS0.cutoff = 0.002
  # matches.samples.keep = read.table("input_data/matches_samples_to_keep.txt", header = F, stringsAsFactors = F)[,1]
  # matches.samples.exclude = read.table("input_data/matches_samples_to_exclude.txt", header = F, stringsAsFactors = F)[,1]
  
  ## Recalculate sample metrics
  leaf_data = leaf_recalculate_sample_metrics(leaf_data = leaf_data)
  
  ## get known relationships
  knownPO = c(apply(leaf_data$sample_data[,c("ID_reps","Parent1")], 1, function(x) paste(sort(x), collapse = "__")), apply(leaf_data$sample_data[,c("ID_reps","Parent2")], 1, function(x) paste(sort(x), collapse = "__")))
  knownPO = grep("Unknown",knownPO, invert = T, value = T)
  knownPO = unique(grep("__",knownPO, value = T))
  
  knownPO2 = knownPO
  ## add Identifier from ID_reps column
  for(pairi in knownPO){
    allIds = sapply(unlist(strsplit(pairi,"__")), function(x) ifelse(!x %in% leaf_data$sample_data$ID, x, list(leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$ID %in% x)])))
    newPairs = unique(c(sapply(c(sapply(allIds[[1]], function(y) paste(y, allIds[[2]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__")),
                        sapply(c(sapply(allIds[[2]], function(y) paste(y, allIds[[1]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__"))
    ))
    knownPO2 = unique(c(knownPO2, newPairs))
    rm(allIds, newPairs)
  }
  knownPO = knownPO2; rm(knownPO2, pairi)
  # //
  
  # sub in IDs
  ibd.robust.coeff = ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.02),] 
  ibd.robust.coeff$ID_reps1 = leaf_data$sample_data$ID_reps[match(ibd.robust.coeff$ID1, leaf_data$sample_data$Sample.Filename)]
  ibd.robust.coeff$ID_reps2 = leaf_data$sample_data$ID_reps[match(ibd.robust.coeff$ID2, leaf_data$sample_data$Sample.Filename)]
  ibd.robust.coeff$pair = ""
  ibd.robust.coeff$pair[which(ibd.robust.coeff$IBS0 < 0.02)] = apply(ibd.robust.coeff[which(ibd.robust.coeff$IBS0 < 0.02),c("ID_reps1","ID_reps2")], 1, function(x) paste(sort(x), collapse = "__"))  
  ibd.robust.coeff$relation = "Unknown"
  ibd.robust.coeff$relation[which(ibd.robust.coeff$pair %in% knownPO)] = "PO"
  
  ## make table with replicate groups
  # dups identified above
  matches.df = ibd.robust.coeff[which(ibd.robust.coeff$kinship > kinship.cutoff & ibd.robust.coeff$IBS0 < IBS0.cutoff),c("ID1","ID2")]
  # sub in IDs
  matches.df$ID1 = leaf_data$sample_data$ID_reps[match(matches.df$ID1, leaf_data$sample_data$Sample.Filename)]
  matches.df$ID2 = leaf_data$sample_data$ID_reps[match(matches.df$ID2, leaf_data$sample_data$Sample.Filename)]
  
  # new table with all reps compiled
  matches.agg.df = data.frame(ID1= unique(c(matches.df[,1], matches.df[,2])), stringsAsFactors = F)
  matches.agg.df$matches = NA
  i=104
  for(i in 1:nrow(matches.agg.df)){
    matches.agg.df$matches[i] = paste(unique(sort(unlist(c(matches.df[which(matches.df$ID1 %in% matches.agg.df$ID1[i] | matches.df$ID2 %in% matches.agg.df$ID1[i]),])))), collapse = "__")
  }
  matches.agg.df = matches.agg.df[which(!duplicated(matches.agg.df$matches)),]
  rm(matches.df)
  
  
  ## process matches.agg.df: sort matches by Call Rate, trim reps, find IDs that are pedigree confirmed by the data
  # 1. sort matches by Call Rate
  matches.agg.df$matchesSorted = unlist(lapply(lapply(sapply(matches.agg.df$matches, function(x) strsplit(x, "__")), function(y) y[order(leaf_data$sample_data$call_rate_recalculated[match(y, leaf_data$sample_data$ID_reps)], decreasing = T)]), function(z) paste(z, collapse="__")))
  
  # 2. trim replicates - if sample is expected replicate then cut from matchesSorted list
  matches.agg.df$matchesSortedTrim = unlist(lapply(lapply(lapply(sapply(matches.agg.df$matchesSorted, function(x) strsplit(x, "__") ), function(y) sapply(y, function(z) unlist(strsplit(z,"_rep"))[1])), function(z) names(z)[which(!duplicated(z))]), function(a) paste(a, collapse="__")))
  
  # 3. find IDs that are pedigree confirmed by the data
  # get table with known and confirmed PO relationships, set looser IBS0 cutoff here to include borderline cases
  ibd.PO.Known.Conf = ibd.robust.coeff[which(ibd.robust.coeff$relation == "PO" & ibd.robust.coeff$IBS0 < 0.006 & ibd.robust.coeff$kinship > 0.1),]
  # ibd.PO.Known.Conf$pair2 = apply(ibd.PO.Known.Conf[,1:2], 1, function(x) paste(sort(gsub("_rep[A-Z]","",x)), collapse = "__"))
  # ibd.PO.Known.Conf = ibd.PO.Known.Conf[which(!duplicated(ibd.PO.Known.Conf$pair2)),]
  
  # Check pedigree
  matches.agg.df$checkPed = FALSE
  matches.agg.df$checkPed[grep("__", matches.agg.df$matchesSortedTrim)] = TRUE
  
  matches.agg.df$IdConfirmed = NA
  matches.agg.df$IdConfirmedNumSupport = NA
  matches.agg.df$checkPedRes = NA
  matches.agg.df$IdConfirmed[which(matches.agg.df$checkPed)] = unlist(lapply(lapply(matches.agg.df$matchesSortedTrim[which(matches.agg.df$checkPed)], function(x) unlist(lapply( strsplit(x, "__"), function(y) sapply(y, function(z) length(grep(z, ibd.PO.Known.Conf$pair))  )))), function(x) ifelse(length(names(x[which(x==max(x) & x != 0)]))>0, names(x[which(x==max(x) & x != 0)]), NA)))
  matches.agg.df$IdConfirmedNumSupport[which(matches.agg.df$checkPed)] = unlist(lapply(lapply(matches.agg.df$matchesSortedTrim[which(matches.agg.df$checkPed)], function(x) unlist(lapply( strsplit(x, "__"), function(y) sapply(y, function(z) length(grep(z, ibd.PO.Known.Conf$pair))  )))), function(x) ifelse(length(names(x[which(x==max(x) & x != 0)]))>0, (x[which(x==max(x) & x != 0)]), NA)))
  matches.agg.df$checkPedRes[which(matches.agg.df$checkPed)] = unlist(lapply(lapply(matches.agg.df$matchesSortedTrim[which(matches.agg.df$checkPed)], function(x) unlist(lapply( strsplit(x, "__"), function(y) sapply(y, function(z) length(grep(z, ibd.PO.Known.Conf$pair))  )))), function(x) paste(x, collapse=",")))
  #
  
  # add numConfirmPOrels to leaf_data$sample_data
  leaf_data$sample_data$numConfirmPOrels = NA
  # numConfirmPOrels = table(gsub("_rep[A-Z]","",c(ibd.PO.Known.Conf$ID1, ibd.PO.Known.Conf$ID2)))
  numConfirmPOrels = table(c(ibd.PO.Known.Conf$ID1, ibd.PO.Known.Conf$ID2))
  leaf_data$sample_data$numConfirmPOrels = numConfirmPOrels[match(leaf_data$sample_data$Sample.Filename, names(numConfirmPOrels))]
  
  # add numFailPOrels to leaf_data$sample_data
  ibd.PO.Known.Fail = ibd.robust.coeff[which(ibd.robust.coeff$relation == "PO" & ibd.robust.coeff$IBS0 > 0.006 ),]
  # ibd.PO.Known.Fail$pair2 = apply(ibd.PO.Known.Fail[,1:2], 1, function(x) paste(sort(gsub("_rep[A-Z]","",x)), collapse = "__"))
  # ibd.PO.Known.Fail = ibd.PO.Known.Fail[which(!duplicated(ibd.PO.Known.Fail$pair2)),]
  leaf_data$sample_data$numFailPOrels = NA
  # numFailPOrels = table(gsub("_rep[A-Z]","",c(ibd.PO.Known.Fail$ID1, ibd.PO.Known.Fail$ID2)))
  numFailPOrels = table(c(ibd.PO.Known.Fail$ID1, ibd.PO.Known.Fail$ID2))
  leaf_data$sample_data$numFailPOrels = numFailPOrels[match(leaf_data$sample_data$Sample.Filename, names(numFailPOrels))]
  #
  
  # Force keep IDs for a set of specified ID_rep's
  matches.agg.df$force.keep = FALSE
  if(!is.null(matches.samples.keep)){
    for(i in 1:length(matches.samples.keep)){
      matches.agg.df$force.keep[grep(matches.samples.keep[i], matches.agg.df$matches)] = TRUE
      matches.agg.df$IdConfirmed[grep(matches.samples.keep[i], matches.agg.df$matches)] = matches.samples.keep[i]
    }
  }
  
  # for cases where pedigree doesn't find correct ID - combine IDs
  matches.agg.df$ambiguous = FALSE
  for(i in 1:nrow(matches.agg.df)){
    if(matches.agg.df$checkPed[i] & is.na(matches.agg.df$IdConfirmed[i])){
      matches.agg.df$IdConfirmed[i] = gsub("__", ";", matches.agg.df$matchesSorted[i])
      matches.agg.df$ambiguous[i] = TRUE
    }
  }
  
  # //
  
  # write.csv(matches.agg.df[,-1],"Sample_Metadata/Identified_Replicates.csv", row.names = F)
  
  # //////////
  
  ############
  ## Add annotation to sampleInfo - info on matches, pedigree hits
  sampleInfo_all_cels = leaf_data$sample_data

    #
  sampleInfo_all_cels$matchOccur = FALSE
  sampleInfo_all_cels$matches = NA
  sampleInfo_all_cels$matchOccurAmbig = NA
  sampleInfo_all_cels$matchIDConfirmed = NA
  sampleInfo_all_cels$ID_reps_amb = sampleInfo_all_cels$ID_reps
  
  #
  i=67
  for(i in 1:nrow(matches.agg.df)){
    hits = unlist(strsplit(matches.agg.df$matchesSorted[i],"__"))
    # hitsTrim = unlist(strsplit(matches.agg.df$matchesSortedTrim[i],"__"))
    sampleInfo_all_cels$matchOccur[which(sampleInfo_all_cels$ID_reps %in% hits)] =  TRUE
    sampleInfo_all_cels$matches[which(sampleInfo_all_cels$ID_reps %in% hits)] = paste(hits, collapse = ";")
    sampleInfo_all_cels$matchOccurAmbig[which(sampleInfo_all_cels$ID_reps %in% hits)] = matches.agg.df$checkPed[i]
    
    if(matches.agg.df$ambiguous[i]){
      sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps == hits[1])] = matches.agg.df$IdConfirmed[i]
      sampleInfo_all_cels$ID_reps_amb[which(sampleInfo_all_cels$ID_reps %in% hits)] = matches.agg.df$IdConfirmed[i]
      sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% hits[-1])] = "ambiguous_ID"
    } else {
      if(matches.agg.df$checkPed[i]){ # cases where pedigree was checked
        if(!is.na(matches.agg.df$IdConfirmed[i])){ # cases where pedigree was confirmed
          sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps == matches.agg.df$IdConfirmed[i])] = matches.agg.df$IdConfirmed[i]
          hits2 = sapply(hits, function(x) unlist(strsplit(x, "_rep"))[1])
          confirm2 = unlist(strsplit(matches.agg.df$IdConfirmed[i], "_rep"))[1]
          falseID = hits[which(!hits2 %in% confirm2)]
          repID = hits[which(!hits %in% c(matches.agg.df$IdConfirmed[i], falseID))]
          sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% falseID)] = "falsified_ID"
          sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% repID)] = "rep_ID"
          rm(hits2, confirm2, falseID, repID)
        } #else { # cases without pedigree info
        # sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% hits)] = "ambiguous_ID"
        # }
      } else { # cases with only replicates, take best call rate
        sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% hits[1])] = hits[1]
        sampleInfo_all_cels$matchIDConfirmed[which(sampleInfo_all_cels$ID_reps %in% hits[-1])] = "rep_ID"
      }
    }
    rm(hits)
  }
  # ///
  
  # Automatically add on suffix for replicates
  sampleInfo_all_cels$rep2 = FALSE
  sampleInfo_all_cels$rep2[which(sampleInfo_all_cels$ID_reps_amb %in% names(table(sampleInfo_all_cels$ID_reps_amb))[which(table(sampleInfo_all_cels$ID_reps_amb) > 1)])] = TRUE
  sampleInfo_all_cels = sampleInfo_all_cels[order(sampleInfo_all_cels$ID_reps_amb),]
  sampleInfo_all_cels$ID_reps_amb[which(sampleInfo_all_cels$rep2)] = unlist(tapply(sampleInfo_all_cels$ID_reps_amb[which(sampleInfo_all_cels$rep2)], sampleInfo_all_cels$ID_reps_amb[which(sampleInfo_all_cels$rep2)], function(x) paste(x,LETTERS[1:length(x)], sep = "_rep")))
  
  
  ## Mark samples to use downstream
  sampleInfo_all_cels$UseDownstream = NA
  sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$matchIDConfirmed %in% c("rep_ID","ambiguous_ID","falsified_ID") | sampleInfo_all_cels$Pass.Fail == "Fail")] = FALSE
  sampleInfo_all_cels$UseDownstream[which((!sampleInfo_all_cels$matchIDConfirmed %in% c("rep_ID","ambiguous_ID","falsified_ID")) & sampleInfo_all_cels$matchOccur & sampleInfo_all_cels$Pass.Fail == "Pass")] = TRUE
  sampleInfo_all_cels$UseDownstream[which(!sampleInfo_all_cels$matchOccur & sampleInfo_all_cels$Pass.Fail == "Pass")] = TRUE
  
  #
  # Exclude
  if(!is.null(matches.samples.exclude)){
    for(i in 1:length(matches.samples.exclude)){
      sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == matches.samples.exclude[i])] = FALSE
    }
  }
  
  # For a few ID, need to keep both replicates, bc each rep was used in GWAS
  if(!is.null(matches.samples.keep)){
    for(i in 1:length(matches.samples.keep)){
      sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == matches.samples.keep[i])] = TRUE
      # print(sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == matches.samples.keep[i])])
    }
  }

  # /////
  
  ## Add columns to sampleInfo
  leaf_data$sample_data = sampleInfo_all_cels[match(colnames(leaf_data$snp_data), rownames(sampleInfo_all_cels)),]
  
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
  
  return(leaf_data)
  
}






#######################################################################
# Pedigree trio analysis
leaf_trio_pedigree = function(leaf_data = NULL, ibd.coeff = NULL, sample.ids = NULL) {
  
  # ibd.coeff = ibd.robust.coeff.filter
  # sample.ids = ids.filter
  # leaf_data = dat4
  
  trio.mendel.errors = list()
  trio.nsnps = list()
  trio.missing.snps = list()
  
  #Loopy loop i=969
  pb <- txtProgressBar(min = 0, max = length(sample.ids), style = 3)
  for(i in 1:length(sample.ids)){
    setTxtProgressBar(pb, i)
    #finds all possible parents of individual i
    poss.par = c(ibd.coeff[which(ibd.coeff$ID1 == sample.ids[i]),2],ibd.coeff[which(ibd.coeff$ID2 == sample.ids[i]),1])
    if(length(poss.par) > 1){ # all(poss.par %in% colnames(snp))
      for(j in 1:(length(poss.par)-1)){
        for(k in (j+1):length(poss.par)){
          # skip cases were Offspring and a Parent are replicates
          # if((!gsub("_rep[A-Z]","", sample.ids[i]) %in% gsub("_rep[A-Z]","", poss.par[j])) & 
          #    (!gsub("_rep[A-Z]","", sample.ids[i]) %in% gsub("_rep[A-Z]","", poss.par[k]))){
            trio = leaf_data$snp_data[,match(c(sample.ids[i],poss.par[j],poss.par[k]), colnames(leaf_data$snp_data))]
            missing.snps = rownames(trio[which(!(trio[,1] != -1 &
                                                   trio[,2] != -1 &
                                                   trio[,3] != -1)
            ),])
            trio = trio[which(trio[,1] != -1 &
                                trio[,2] != -1 &
                                trio[,3] != -1
            ),]
            mendel.errors = rownames(trio[which(
              (trio[,1] == 1 & trio[,2] == 0 & trio[,3] == 0) |
                (trio[,1] == 1 & trio[,2] == 2 & trio[,3] == 2) |
                (trio[,1] == 0 & (trio[,2] == 2 | trio[,3] == 2)) |
                (trio[,1] == 2 & (trio[,2] == 0 | trio[,3] == 0))  
            ),])
            
            trio.mendel.errors = c(trio.mendel.errors, list(mendel.errors))
            names(trio.mendel.errors)[length(trio.mendel.errors)] = paste(colnames(trio), collapse ="__")
            
            trio.nsnps = c(trio.nsnps, list(nrow(trio)))
            names(trio.nsnps)[length(trio.nsnps)] = paste(colnames(trio), collapse ="__")
            
            trio.missing.snps = c(trio.missing.snps, list(missing.snps))
            names(trio.missing.snps)[length(trio.missing.snps)] = paste(colnames(trio), collapse ="__")
            
            # print(c(i,j,k))
            rm(mendel.errors, missing.snps, trio)
          # }
        } #end k loop
      } # end j loop
    }
  } # end i loop
  close(pb)
  # //
  
  # trio-wise
  trio.transgressions = unlist(lapply(trio.mendel.errors, length))
  #
  trio.df = data.frame(trio=names(trio.transgressions), transgressions=trio.transgressions, total.tests = unlist(trio.nsnps),  stringsAsFactors = F)
  trio.df$ID = sapply(trio.df$trio, function(x) unlist(strsplit(x, split = "__"))[1])
  trio.df$Parent1 = sapply(trio.df$trio, function(x) unlist(strsplit(x, split = "__"))[2])
  trio.df$Parent2 = sapply(trio.df$trio, function(x) unlist(strsplit(x, split = "__"))[3])
  trio.df$ID.name = leaf_data$sample_data$ID_reps_amb[match(trio.df$ID, leaf_data$sample_data$Sample.Filename)]
  trio.df$Parent1.name = leaf_data$sample_data$ID_reps_amb[match(trio.df$Parent1, leaf_data$sample_data$Sample.Filename)]
  trio.df$Parent2.name = leaf_data$sample_data$ID_reps_amb[match(trio.df$Parent2, leaf_data$sample_data$Sample.Filename)]

  #
  trio.df$ttr = trio.df$transgressions / trio.df$total.tests
  #
  
  return(trio.df)
}







###########################################################
# Process IBD-KING results
leaf_duo_ibd_king_process = function(leaf_data = NULL, ibd.coeff = NULL, IBS0.cutoff.max = 0.002) {
  
  # ibd.coeff = ibd.robust.coeff.filter
  # leaf_data = dat1
  
  ##
  sampleInfo.ped = leaf_data$sample_data[which(leaf_data$sample_data$UseDownstream),c("ID","ID_reps","ID_reps_amb","Parent1","Parent2")]
  sampleInfo.ped$parentsPair = apply(sampleInfo.ped[,c("Parent1","Parent2")], 1, function(x) paste(sort(x), collapse = "__"))
  # get replicates
  replicatePairs = grep("__", 
                        tapply(sampleInfo.ped$ID_reps_amb, sampleInfo.ped$ID, function(x) paste(sort(x), collapse = "__")), 
                        value=T)
  
  # full-sibs
  familiesList = tapply(sampleInfo.ped$ID_reps_amb, sampleInfo.ped$parentsPair, c)
  familiesList = familiesList[which(names(familiesList) != "")]
  familiesList = familiesList[which(names(familiesList) != "UnknownA__UnknownB")]
  familiesList = familiesList[which(lapply(familiesList, length) > 1)]
  familiesList = lapply(familiesList, function(x) { xx = data.frame(t(combn(x, 2)), stringsAsFactors = F); 
  apply(xx, 1, function(y) paste(sort(y), collapse = "__"))})
  fullsibs = unlist(familiesList)
  # length(fullsibs)
  fullsibs = setdiff(fullsibs, replicatePairs) # remove replicatePairs
  
  #
  # half-sibs
  ped1 = sampleInfo.ped[,c("ID_reps_amb","Parent1")]
  ped2 = sampleInfo.ped[,c("ID_reps_amb","Parent2")]
  colnames(ped2)[2] = "Parent1"
  ped1 = rbind(ped1, ped2); rm(ped2)
  poList = tapply(ped1$ID, ped1$Parent1, c)
  poList = poList[which(lapply(poList, length) > 1)]
  poList = poList[-grep("Unknown", names(poList))]
  poList = lapply(poList, function(x) { xx = data.frame(t(combn(x, 2)), stringsAsFactors = F); 
  apply(xx, 1, function(y) paste(sort(y), collapse = "__"))})
  poList = unlist(poList)
  halfsibs = setdiff(poList, fullsibs)
  halfsibs = setdiff(halfsibs, replicatePairs)
  
  # POs
  poKnown = apply(ped1, 1, function(y) paste(sort(y), collapse = "__"))
  poKnown = grep("Unknown",  poKnown, value=T, invert = T)
  poKnown = grep("__",  poKnown, value=T)
  rm(ped1, familiesList)
  
  poKnown2 = poKnown
  ## add Identifier from ID_reps column
  for(pairi in poKnown){
    allIds = sapply(unlist(strsplit(pairi,"__")), function(x) ifelse(!x %in% leaf_data$sample_data$ID, x, list(leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$ID %in% x)])))
    newPairs = unique(c(sapply(c(sapply(allIds[[1]], function(y) paste(y, allIds[[2]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__")),
                        sapply(c(sapply(allIds[[2]], function(y) paste(y, allIds[[1]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__"))
    ))
    poKnown2 = unique(c(poKnown2, newPairs))
    rm(allIds, newPairs)
  }
  # tail(poKnown2, 20)
  poKnown = poKnown2; rm(poKnown2, pairi)
  
  
  ##
  ibd.coeff$pair = apply(ibd.coeff[,c(1,2)], 1, function(x) paste(sort(
    leaf_data$sample_data$ID_reps_amb[match(x, leaf_data$sample_data$Sample.Filename)]
  ), collapse = "__"))
  
  ibd.coeff$relation = "Unknown"
  ibd.coeff$relation[which(ibd.coeff$pair %in% replicatePairs)] = "Replicate pair"
  ibd.coeff$relation[which(ibd.coeff$pair %in% poKnown)] = "PO - recorded"
  ibd.coeff$relation[which(ibd.coeff$pair %in% fullsibs)] = "Fullsibs - recorded"
  ibd.coeff$relation[which(ibd.coeff$pair %in% halfsibs)] = "Halfsibs - recorded"
  
  
  ########
  # Make calls (SNP error rate from technical replicates = 0.00099 and 0.00107)
  ibd.coeff$KingCall = "Unknown"
  ibd.coeff$KingCall[which(ibd.coeff$IBS0 < IBS0.cutoff.max & 
                             ibd.coeff$kinship > 0.4 &
                             ibd.coeff$relation == "replicate pair")] = "Replicate pair - Confirmed"
  ibd.coeff$KingCall[which(ibd.coeff$IBS0 < IBS0.cutoff.max & 
                             # ibd.coeff$kinship > 0.1 & 
                             ibd.coeff$kinship < 0.4 & 
                             ibd.coeff$relation != "PO - recorded")] = "PO called - New"
  ibd.coeff$KingCall[which(ibd.coeff$IBS0 < IBS0.cutoff.max & 
                             # ibd.coeff$kinship > 0.1 &
                             ibd.coeff$kinship < 0.4 & 
                             ibd.coeff$relation == "PO - recorded")] = "PO called - Confirmed"
  ibd.coeff$KingCall[which(!(ibd.coeff$IBS0 < IBS0.cutoff.max 
                             # & ibd.coeff$kinship > 0.1
                             ) & ibd.coeff$relation == "PO - recorded")] = "PO failed"
  ibd.coeff$KingCall[which((ibd.coeff$IBS0 < 0.01 
                            # & ibd.coeff$kinship > 0.1
                            ) & ibd.coeff$KingCall == "PO failed")] = "PO failed - borderline"
  ibd.coeff$KingCall[which(ibd.coeff$relation == "Fullsibs - recorded")] = "Fullsibs - recorded"
  ibd.coeff$KingCall[which(ibd.coeff$relation == "Halfsibs - recorded")] = "Halfsibs - recorded"
  
  return(ibd.coeff)
}



######################################################################
# Get pedigree from IBD-KING
leaf_duo_ibd_king_pedigree = function(leaf_data = NULL, ibd.coeff = NULL, kinship.cutoff.min = 0.1, kinship.cutoff.max = 0.4, IBS0.cutoff.max = 0.002) {
  
  # leaf_data = dat1
  # ibd.coeff = ibd.robust.coeff.filter.res
  # kinship.cutoff.min = 0.1
  # kinship.cutoff.max = 0.4
  # IBS0.cutoff.max = 0.002
  
  ## Get putative PO pairs
  putPO = ibd.coeff[which((ibd.coeff$IBS0 < IBS0.cutoff.max & 
                             ibd.coeff$kinship > kinship.cutoff.min & 
                             ibd.coeff$kinship < kinship.cutoff.max) | 
                            ibd.coeff$KingCall == "PO failed - borderline"),]
  ## Parse putPO object
  putPO$pair.filename = apply(putPO[,c("ID1","ID2")], 1, function(x) paste(sort(x), collapse = "__"))
  # ifelse to catch cases without Year metadata
  putPO$Offspring = apply(putPO[,c(1,2)], 1, function(x) {
    ifelse(
      length(which.max(leaf_data$sample_data$Year[match(x,leaf_data$sample_data$Sample.Filename)])) > 0, x[which.max(leaf_data$sample_data$Year[match(x,leaf_data$sample_data$Sample.Filename)])], x[1])
  }
  )
  putPO$Parent = apply(putPO[,c(1,2)], 1, function(x) {
    ifelse(
      length(which.min(leaf_data$sample_data$Year[match(x,leaf_data$sample_data$Sample.Filename)])) > 0, x[which.min(leaf_data$sample_data$Year[match(x,leaf_data$sample_data$Sample.Filename)])], x[2])
  }
  )
  putPO$flag = "Pass"
  putPO$flag[which(apply(putPO[,c(1,2)], 1, function(x) length(unique(leaf_data$sample_data$Year[match(x,leaf_data$sample_data$Sample.Filename)])) == 1))] = "SameYear"
  putPO$Offspring[which(putPO$flag == "SameYear")] = putPO$ID1[which(putPO$flag == "SameYear")]
  putPO$Parent[which(putPO$flag == "SameYear")] = putPO$ID2[which(putPO$flag == "SameYear")]

  putPed = data.frame(ID = unique(putPO$Offspring), stringsAsFactors = F)
  putPed$Parents = sapply(putPed$ID, function(x) paste(putPO$Parent[which(putPO$Offspring == x)],collapse = "__"))
  putPed$numParents = sapply(putPed$ID, function(x) length(putPO$Parent[which(putPO$Offspring == x)]))
  putPed$Parents.Rel = sapply(putPed$ID, function(x) paste(putPO$relation[which(putPO$Offspring == x)],collapse = "__"))
  putPed$Parents.ibs0 = sapply(putPed$ID, function(x) paste(round(putPO$IBS0[which(putPO$Offspring == x)],10),collapse = "__"))
  putPed$Parents.kin = sapply(putPed$ID, function(x) paste(round(putPO$kinship[which(putPO$Offspring == x)],4),collapse = "__"))
  
  putPed$Parent1 = sapply(putPed$Parents, function(x) unlist(strsplit(x, split="__"))[1])
  putPed$Parent2 = sapply(putPed$Parents, function(x) unlist(strsplit(x, split="__"))[2])
  putPed$Parent1.IBS0 = sapply(putPed$Parents.ibs0, function(x) unlist(strsplit(x, split="__"))[1])
  putPed$Parent2.IBS0 = sapply(putPed$Parents.ibs0, function(x) unlist(strsplit(x, split="__"))[2])
  # rownames(putPed) = putPed$ID
  
  putPed$ID.name = leaf_data$sample_data$ID_reps_amb[match(putPed$ID, leaf_data$sample_data$Sample.Filename)]
  putPed$Parent1.name = leaf_data$sample_data$ID_reps_amb[match(putPed$Parent1, leaf_data$sample_data$Sample.Filename)]
  putPed$Parent2.name = leaf_data$sample_data$ID_reps_amb[match(putPed$Parent2, leaf_data$sample_data$Sample.Filename)]
  
  putPed$flag = apply(putPed[,c("ID","Parent1","Parent2")], 1, function(x){
    ifelse(paste(sort(x[1:2]), collapse = "__") %in% putPO$pair.filename[which(putPO$flag == "SameYear")], "SameYear", NA)
  })
    
  
  results = list(duo.results.pairs = putPO, duo.results.pedigree = putPed)
  return(results)
  
}



######################################################################
# Consolidate two pedigree dataframes
leaf_consolidate_pedigrees = function(ped1 = NULL, ped2 = NULL, resolve.conflict = NULL, poFailed = NULL) {
  # rename colnames
  colnames(ped1) = c("ID","Parent1","Parent2")
  colnames(ped2) = c("ID","Parent1","Parent2")
  
  # new consolidated pedigree
  ped = data.frame(ID = unique(c(ped1$ID, ped2$ID)), con.Parent1 = NA, con.Parent2 = NA, ped1.Parent1 = NA, ped1.Parent2 = NA, ped2.Parent1 = NA, ped2.Parent2 = NA, stringsAsFactors = F)
  ped$ped1.Parent1 = ped1$Parent1[match(ped$ID, ped1$ID)]
  ped$ped1.Parent2 = ped1$Parent2[match(ped$ID, ped1$ID)]
  ped$ped2.Parent1 = ped2$Parent1[match(ped$ID, ped2$ID)]
  ped$ped2.Parent2 = ped2$Parent2[match(ped$ID, ped2$ID)]
  
  # clean up cases where parent1 empty but parent2 not empty
  for(i in 1:nrow(ped)){
    if(is.na(ped$ped1.Parent1[i]) & (!is.na(ped$ped1.Parent2[i])) ){
      hold = ped$ped1.Parent2[i]
      ped$ped1.Parent1[i] = hold
      ped$ped1.Parent2[i] = NA
    }
    if(is.na(ped$ped2.Parent1[i]) & (!is.na(ped$ped2.Parent2[i])) ){
      hold = ped$ped2.Parent2[i]
      ped$ped2.Parent1[i] = hold
      ped$ped2.Parent2[i] = NA
    }
  }
  
  ped$notes = NA
  
  #
  for(i in 1:nrow(ped)){
    pedi = as.character(ped[i,c("ped1.Parent1","ped1.Parent2","ped2.Parent1","ped2.Parent2")])

    # 2/2 matches
    if((length(na.omit(pedi[1:2])) == 2 & length(na.omit(pedi[3:4])) == 2) &
       identical(sort(pedi[1:2]), sort(pedi[3:4])) )
      {
      ped$con.Parent1[i] = pedi[1]
      ped$con.Parent2[i] = pedi[2]
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "2/2 matches")), collapse = ";")
    } else if((length(na.omit(pedi[1:2])) == 2 & length(na.omit(pedi[3:4])) == 2) & 
              (! identical(sort(pedi[1:2]), sort(pedi[3:4])) ) ) {
      if(length(table(pedi)) == 4){
        ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/2 match - ped1 CONFLICT with ped2")), collapse = ";")
      }
      if(length(table(pedi)) == 3){
        ped$notes[i] = paste(na.omit(c(ped$notes[i], "1/2 match - ped1 CONFLICT with ped2")), collapse = ";")
      }
      if(resolve.conflict == "ped1"){
        ped$con.Parent1[i] = pedi[1]
        ped$con.Parent2[i] = pedi[2]
      }
      if(resolve.conflict == "ped2"){
        ped$con.Parent1[i] = pedi[3]
        ped$con.Parent2[i] = pedi[4]
      }
    }
    
    # 1/1 match - ped2 has only 1 parent
    if((length(na.omit(pedi[1:2])) == 2 & length(na.omit(pedi[3:4])) == 1) & (
       pedi[1] %in% pedi[3:4] |
       pedi[2] %in% pedi[3:4] )
    ){
      ped$con.Parent1[i] = pedi[1]
      ped$con.Parent2[i] = pedi[2]
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "1/1 match - ped2 has only 1 parent")), collapse = ";")
    } else if((length(na.omit(pedi[1:2])) == 2 & length(na.omit(pedi[3:4])) == 1) & (!(
      pedi[1] %in% pedi[3:4] |
      pedi[2] %in% pedi[3:4]) )
    ){
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/1 match - ped2 CONFLICT with ped1")), collapse = ";")
      if(resolve.conflict == "ped1"){
        ped$con.Parent1[i] = pedi[1]
        ped$con.Parent2[i] = pedi[2]
      }
      if(resolve.conflict == "ped2"){
        ped$con.Parent1[i] = pedi[3]
        ped$con.Parent2[i] = pedi[4]
      }
    }
    
    # 1/1 match - ped1 has only 1 parent
    if((length(na.omit(pedi[1:2])) == 1 & length(na.omit(pedi[3:4])) == 2) & (
      pedi[3] %in% pedi[1:2] |
      pedi[4] %in% pedi[1:2] )
    ){
      ped$con.Parent1[i] = pedi[3]
      ped$con.Parent2[i] = pedi[4]
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "1/1 match - ped1 has only 1 parent")), collapse = ";")
    } else if((length(na.omit(pedi[1:2])) == 1 & length(na.omit(pedi[3:4])) == 2) & (!(
      pedi[3] %in% pedi[1:2] |
      pedi[4] %in% pedi[1:2] ))
    ){
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/1 match - ped1 CONFLICT with ped2")), collapse = ";")
      if(resolve.conflict == "ped1"){
        ped$con.Parent1[i] = pedi[1]
        ped$con.Parent2[i] = pedi[2]
      }
      if(resolve.conflict == "ped2"){
        ped$con.Parent1[i] = pedi[3]
        ped$con.Parent2[i] = pedi[4]
      }
    }
    
    # 0/0 matches - ped1 has 0 parents
    if(length(na.omit(pedi[1:2])) == 0){
      ped$con.Parent1[i] = pedi[3]
      ped$con.Parent2[i] = pedi[4]
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/0 matches - ped1 has 0 parents")), collapse = ";")
    }
    
    # 0/0 matches - ped2 has 0 parents
    if(length(na.omit(pedi[3:4])) == 0){
      ped$con.Parent1[i] = pedi[1]
      ped$con.Parent2[i] = pedi[2]
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/0 matches - ped2 has 0 parents")), collapse = ";")
    }
    
    # 0/1 matches - CONFLICT - ped1 and ped2 each have 1 parent
    if(length(na.omit(pedi[1:2])) == 1 & 
       length(na.omit(pedi[3:4])) == 1 & 
       pedi[1] != pedi[3]){
      ped$notes[i] = paste(na.omit(c(ped$notes[i], "0/1 matches - CONFLICT - ped1 and ped2 each have 1 parent")), collapse = ";")
      if(resolve.conflict == "ped1"){
        ped$con.Parent1[i] = pedi[1]
        ped$con.Parent2[i] = pedi[2]
      }
      if(resolve.conflict == "ped2"){
        ped$con.Parent1[i] = pedi[3]
        ped$con.Parent2[i] = pedi[4]
      }
    }
    
    
  } # end of for loop
  #
  # ped = ped[-which(ped$notes == "2/2 matches"),]
  if(!is.null(poFailed)){
    for(i in 1:nrow(ped)){
      if(paste(sort(c(ped$ID[i],ped$con.Parent1[i])), collapse = "__") %in% poFailed)
      {
        ped$con.Parent1[i] = NA
        ped$notes[i] = paste(na.omit(c(ped$notes[i], "parent1 FAIL by KING")), collapse = ";")
      }
      if(paste(sort(c(ped$ID[i],ped$con.Parent2[i])), collapse = "__") %in% poFailed)
      {
        ped$con.Parent2[i] = NA
        ped$notes[i] = paste(na.omit(c(ped$notes[i], "parent2 FAIL by KING")), collapse = ";")
      }
    }
  }
  
  return(ped)
}

#######################################################