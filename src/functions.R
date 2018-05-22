#Functions for the Leaf_print project
#Tom Poorten, Mitchell Feldmann, Randi Famula

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
  
  # leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$rep)] = unlist(tapply(leaf_data$sample_data$ID[which(leaf_data$sample_data$rep)], leaf_data$sample_data$ID[which(leaf_data$sample_data$rep)], function(x) paste(x,LETTERS[1:length(x)], sep = "_rep")))
  
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
leaf_analyze_matches = function(leaf_data = NULL, ibd.robust.coeff = NULL, kinship.cutoff = 0.44, IBS0.cutoff = 0.002){
  # leaf_data = dat3
  # kinship.cutoff = 0.44
  # IBS0.cutoff = 0.002
  
  ## get known relationships
  knownPO = c(apply(leaf_data$sample_data[,c("ID_reps","Parent1")], 1, function(x) paste(sort(x), collapse = "__")), apply(leaf_data$sample_data[,c("ID_reps","Parent2")], 1, function(x) paste(sort(x), collapse = "__")))
  knownPO = grep("Unknown",knownPO, invert = T, value = T)
  knownPO = unique(grep("__",knownPO, value = T))
  
  knownPO2 = knownPO
  # pairi = "63C125P039__75C034P105"
  # pairi = "60C019P002__Rockhill_Derivative"
  ## add Identifier from ID_reps column
  for(pairi in knownPO){
    allIds = sapply(unlist(strsplit(pairi,"__")), function(x) ifelse(!x %in% leaf_data$sample_data$ID, x, list(leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$ID %in% x)])))
    newPairs = unique(c(sapply(c(sapply(allIds[[1]], function(y) paste(y, allIds[[2]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__")),
                        sapply(c(sapply(allIds[[2]], function(y) paste(y, allIds[[1]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__"))
    ))
    knownPO2 = unique(c(knownPO2, newPairs))
    rm(allIds, newPairs)
  }
  # tail(knownPO2, 20)
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
  # update Tribute - USDA ped = EB18 x MDUS 4258; emperical ped = MDUS_5189 x Everbearing_372
  # PI551863-Himiko - previous analysis, I removed this sample with gsub, not sure why
  # get table with known and confirmed PO relationships, set looser IBS0 cutoff here to include borderline cases
  ibd.PO.Known.Conf = ibd.robust.coeff[which(ibd.robust.coeff$relation == "PO" & ibd.robust.coeff$IBS0 < 0.006 & ibd.robust.coeff$kinship > 0.1),]
  ibd.PO.Known.Conf$pair2 = apply(ibd.PO.Known.Conf[,1:2], 1, function(x) paste(sort(gsub("_rep[A-Z]","",x)), collapse = "__"))
  ibd.PO.Known.Conf = ibd.PO.Known.Conf[which(!duplicated(ibd.PO.Known.Conf$pair2)),]
  
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
  numConfirmPOrels = table(gsub("_rep[A-Z]","",c(ibd.PO.Known.Conf$ID1, ibd.PO.Known.Conf$ID2)))
  leaf_data$sample_data$numConfirmPOrels = numConfirmPOrels[match(leaf_data$sample_data$Sample.Filename, names(numConfirmPOrels))]
  # add numFailPOrels to leaf_data$sample_data
  ibd.PO.Known.Fail = ibd.robust.coeff[which(ibd.robust.coeff$relation == "PO" & ibd.robust.coeff$IBS0 > 0.006 ),]
  ibd.PO.Known.Fail$pair2 = apply(ibd.PO.Known.Fail[,1:2], 1, function(x) paste(sort(gsub("_rep[A-Z]","",x)), collapse = "__"))
  ibd.PO.Known.Fail = ibd.PO.Known.Fail[which(!duplicated(ibd.PO.Known.Fail$pair2)),]
  leaf_data$sample_data$numFailPOrels = NA
  numFailPOrels = table(gsub("_rep[A-Z]","",c(ibd.PO.Known.Fail$ID1, ibd.PO.Known.Fail$ID2)))
  leaf_data$sample_data$numFailPOrels = numFailPOrels[match(leaf_data$sample_data$Sample.Filename, names(numFailPOrels))]
  #
  
  # Manual edits
  # 1. I am pretty sure what is called 05C137P006 is actually 39C082P019. I think that was part of the pipetting mistake from Batch 3. 
  #       I don’t know if there is any way for you to check it or if it just should be dropped from the dataset.
  # 2. 
  # I went through this list. I clearly made a mistake in the submission of:
  # 05C137P006- not on Tom’s list, but I am almost positive that I submitted 10C157P001 instead
  # 58C046P002_1AAAB (really submitted Monterey) # AUTOMATICALLY CORRECTED BY PEDIGREE CHECKING
  # PI551529_1BBAA (really submitted San Andreas) # AUTOMATICALLY CORRECTED BY PEDIGREE CHECKING
  # PI551614_1BBAA (really submitted Portola) # AUTOMATICALLY CORRECTED BY PEDIGREE CHECKING
  # PI551799_1BBAA (really submitted Cabrillo) # AUTOMATICALLY CORRECTED BY PEDIGREE CHECKING
  # make manual edits
  matches.agg.df$IdConfirmed[grep("05C137P006_repB", matches.agg.df$matches)] = "39C082P019"
  matches.agg.df$IdConfirmed[grep("37C020P045_repB", matches.agg.df$matches)] = "37C020P045_repB"
  matches.agg.df$IdConfirmed[grep("64C028P018", matches.agg.df$matches)] = "64C028P018"
  matches.agg.df$IdConfirmed[grep("Sitka", matches.agg.df$matches)] = "Sitka"
  matches.agg.df$IdConfirmed[grep("Tribute", matches.agg.df$matches)] = "Tribute"
  matches.agg.df$IdConfirmed[grep("Sparkle", matches.agg.df$matches)] = "Sparkle"
  matches.agg.df$IdConfirmed[grep("Tago", matches.agg.df$matches)] = "Tago"
  matches.agg.df$IdConfirmed[grep("Gorella", matches.agg.df$matches)] = "Gorella"
  matches.agg.df$IdConfirmed[grep("11C018P001", matches.agg.df$matches)] = "11C018P001" # random pick of 1 of 2 fullsibs?
  matches.agg.df$IdConfirmed[grep("Midway", matches.agg.df$matches)] = "Midway" # matches Early_Midway, so keep 1 -> Midway
  
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
  # sampleInfo_all_cels = read.csv("Sample_Metadata/sampleInfo_all_cels.csv", stringsAsFactors = F, header=T)
  sampleInfo_all_cels = leaf_data$sample_data
  #
  # colnames(sampleInfo_all_cels)
  # sampleInfo_all_cels2 = sampleInfo_all_cels[,c(9:14,39)]
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
  
  # Get only table of known replicates
  sampleInfo_reps = sampleInfo_all_cels[which(sampleInfo_all_cels$rep),] 
  #
  # Check by eye
  # Brighton - 72C271P105_repA, 72C271P105_repB - do not match
  #       - 72C271P105_repA is pedigree confirmed as Brighton, although IBD numbers are borderline
  #       - 72C271P105_repB is falsified as Brighton
  # Edit manually
  # 
  sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == "72C271P105_repA")] = TRUE
  sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == "72C271P105_repB")] = FALSE
  
  # For a few ID, need to keep both replicates, bc each rep was used in GWAS
  sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == "51S001P001_repB")] = TRUE
  sampleInfo_all_cels$UseDownstream[which(sampleInfo_all_cels$ID_reps == "65CA1501_repA")] = TRUE
  
  # /////
  
  ## Add columns to sampleInfo
  # colnames(sampleInfo_all_cels)
  # sampleInfo = cbind(sampleInfo, sampleInfo_all_cels[match(sampleInfo$Sample.Filename, sampleInfo_all_cels$Sample.Filename), 
                                                     # c("matchOccur","matches","matchOccurAmbig","matchIDConfirmed","ID_reps_amb","UseDownstream")])
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
  for(i in 1:length(sample.ids)){
    
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
            
            print(c(i,j,k))
            rm(mendel.errors, missing.snps, trio)
          # }
        } #end k loop
      } # end j loop
    }
  } # end i loop
  
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
leaf_duo_ibd_king_process = function(leaf_data = NULL, ibd.coeff = NULL) {
  
  # ibd.coeff = ibd.robust.coeff.filter
  # leaf_data = dat4
  
  ##
  sampleInfo.ped = leaf_data$sample_data[which(leaf_data$sample_data$UseDownstream),c("ID","ID_reps","ID_reps_amb","Parent1","Parent2")]
  sampleInfo.ped$parentsPair = apply(sampleInfo.ped[,c("Parent1","Parent2")], 1, function(x) paste(sort(x), collapse = "__"))
  # full-sibs
  familiesList = tapply(sampleInfo.ped$ID_reps_amb, sampleInfo.ped$parentsPair, c)
  familiesList = familiesList[which(names(familiesList) != "")]
  familiesList = familiesList[which(names(familiesList) != "UnknownA__UnknownB")]
  familiesList = familiesList[which(lapply(familiesList, length) > 1)]
  familiesList = lapply(familiesList, function(x) { xx = data.frame(t(combn(x, 2)), stringsAsFactors = F); 
  apply(xx, 1, function(y) paste(sort(y), collapse = "__"))})
  fullsibs = unlist(familiesList)
  length(fullsibs)
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
  halfsibs = poList[which(!poList %in% fullsibs)]
  
  # POs
  poCalled = apply(ped1, 1, function(y) paste(sort(y), collapse = "__"))
  poCalled = grep("Unknown",  poCalled, value=T, invert = T)
  poCalled = grep("__",  poCalled, value=T)
  rm(ped1, familiesList)
  
  poCalled2 = poCalled
  ## add Identifier from ID_reps column
  for(pairi in poCalled){
    allIds = sapply(unlist(strsplit(pairi,"__")), function(x) ifelse(!x %in% leaf_data$sample_data$ID, x, list(leaf_data$sample_data$ID_reps[which(leaf_data$sample_data$ID %in% x)])))
    newPairs = unique(c(sapply(c(sapply(allIds[[1]], function(y) paste(y, allIds[[2]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__")),
                        sapply(c(sapply(allIds[[2]], function(y) paste(y, allIds[[1]], sep="__"))), function(z) paste(sort(c(unlist(strsplit(z, "__")))),collapse="__"))
    ))
    poCalled2 = unique(c(poCalled2, newPairs))
    rm(allIds, newPairs)
  }
  # tail(poCalled2, 20)
  poCalled = poCalled2; rm(poCalled2, pairi)
  
  
  
  ##
  ibd.coeff$pair = apply(ibd.coeff[,c(1,2)], 1, function(x) paste(sort(
    leaf_data$sample_data$ID_reps_amb[match(x, leaf_data$sample_data$Sample.Filename)]
  ), collapse = "__"))
  
  ibd.coeff$relation = "Unknown"
  ibd.coeff$relation[which(ibd.coeff$pair %in% poCalled)] = "PO - recorded"
  ibd.coeff$relation[which(ibd.coeff$pair %in% fullsibs)] = "Fullsibs - recorded"
  ibd.coeff$relation[which(ibd.coeff$pair %in% halfsibs)] = "Halfsibs - recorded"
  
  
  ########
  # Make calls (SNP error rate from technical replicates = 0.00099 and 0.00107)
  ibd.coeff$KingCall = "Unknown"
  ibd.coeff$KingCall[which(ibd.coeff$IBS0 < 0.002 & ibd.coeff$kinship > 0.1 & ibd.coeff$relation != "PO - recorded")] = "PO called - New"
  ibd.coeff$KingCall[which(ibd.coeff$IBS0 < 0.002 & ibd.coeff$kinship > 0.1 & ibd.coeff$relation == "PO - recorded")] = "PO called - Confirmed"
  ibd.coeff$KingCall[which(!(ibd.coeff$IBS0 < 0.002 & ibd.coeff$kinship > 0.1) & ibd.coeff$relation == "PO - recorded")] = "PO failed"
  ibd.coeff$KingCall[which((ibd.coeff$IBS0 < 0.01 & ibd.coeff$kinship > 0.1) & ibd.coeff$KingCall == "PO failed")] = "PO failed - borderline"
  ibd.coeff$KingCall[which(ibd.coeff$relation == "Fullsibs - recorded")] = "Fullsibs - recorded"
  ibd.coeff$KingCall[which(ibd.coeff$relation == "Halfsibs - recorded")] = "Halfsibs - recorded"
  
  return(ibd.coeff)
}



######################################################################
# Get pedigree from IBD-KING
leaf_duo_ibd_king_pedigree = function(leaf_data = NULL, ibd.coeff = NULL, kinship.cutoff.min = 0.1, kinship.cutoff.max = 0.4, IBS0.cutoff.max = 0.002) {
  
  # leaf_data = dat4
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
  # putPOall = putPO
  # putPO = putPO[which(putPO$flag == "Pass"),]
  # putPO$summary = paste0("IBS0:",round(putPO$IBS0,5),",kinship:",round(putPO$kinship,5))
  
  putPed = data.frame(ID = unique(putPO$Offspring), stringsAsFactors = F)
  putPed$Parents = sapply(putPed$ID, function(x) paste(putPO$Parent[which(putPO$Offspring == x)],collapse = "__"))
  putPed$numParents = sapply(putPed$ID, function(x) length(putPO$Parent[which(putPO$Offspring == x)]))
  putPed$Parents.Rel = sapply(putPed$ID, function(x) paste(putPO$relation[which(putPO$Offspring == x)],collapse = "__"))
  putPed$Parents.ibs0 = sapply(putPed$ID, function(x) paste(round(putPO$IBS0[which(putPO$Offspring == x)],10),collapse = "__"))
  putPed$Parents.kin = sapply(putPed$ID, function(x) paste(round(putPO$kinship[which(putPO$Offspring == x)],4),collapse = "__"))
  
  # putPed2 = putPed[which(putPed$numParents <= 2),]
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
  
  # write.csv(putPO, "Pedigree/putative_pedigree_KING_duo_results.csv", row.names = F)
  # write.csv(putPed, "Pedigree/putative_pedigree_KING_trioConversion_results.csv", row.names = F)
  # 
  # putPO.borderline = ibd.coeff[which((ibd.coeff$IBS0 < 0.01 & ibd.coeff$kinship > .1 & ibd.coeff$kinship < .4) | ibd.coeff$KingCall == "PO failed - borderline"),]
  # write.csv(putPO.borderline, "Pedigree/KING_results_cut0.01.csv", row.names = F)
}



