#' @export
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")



#' @export
rainfall_conv = c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G", "InDel")
names(rainfall_conv) = c('A>G', 'T>C', 'C>T', 'G>A', 'A>T', 'T>A', 'A>C', 'T>G', 'C>A', 'G>T', 'C>G', 'G>C', 'InDel')



ssh_session <<- NULL



#' @export
cache_output = function(result_df,
                        function_name,
                        clobber_mode = F,
                        get_existing = F,
                        function_params = list(region = "chr3:98300000-198022430", bin_size=2000, seq_type="genome"),
                        additional_details = list(foreground = "DLBCL_FL_BL", background = "CLL_MM_MCL")){

  cache_file_name = paste0(check_config_value(config::get("repo_base")),"cached_results/", function_name)
  for (param in names(function_params)[order(names(function_params))]){
    cache_file_name = paste0(cache_file_name,"--",param,"-",function_params[[param]])
  }
  for (detail in names(additional_details)){
    cache_file_name = paste0(cache_file_name,"--",detail,"-",additional_details[[detail]])
  }
  cache_file_name = paste0(cache_file_name,".tsv")
  if(file.exists(cache_file_name)){
    if(get_existing){
      result_df = suppressMessages(read_tsv(cache_file_name))
      return(result_df)
    }
    if(!clobber_mode){
      warning(paste("file",cache_file_name,"exists!"))
      stop("Will not overwrite unless you rerun this in clobber_mode = TRUE")
    }
  }else{
    if(get_existing){
      stop(paste("cannot find cached result for this parameter combination",cache_file_name))
    }
  }

  message(paste("creating/overwriting",cache_file_name))
  write_tsv(result_df,file=cache_file_name)
}



#' @title Collate Ancestry.
#'
#' @description Gather ancestry information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param somalier_output Somalier ancestery.tsv
#'
#' @return A table.
#'
#' @import stringr readr dplyr
#'
#' @noRd
#'
#' @examples
#' table = collate_ancestry(sample_table = "my_sample_table.txt")
#'
#' @export
collate_ancestry = function(sample_table,
                            seq_type_filter="genome",
                            somalier_output){
  if(seq_type_filter=="capture"){
    message("skipping ancestry for this seq_type")
    return(sample_table)
  }
  if(missing(somalier_output)){
    somalier_output = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/results/gambl/somalier_current/02-ancestry/2020_08_07.somalier-ancestry.tsv"
  }
  somalier_all = suppressMessages(read_tsv(somalier_output))
  somalier_all = mutate(somalier_all, sample_id = str_remove(`#sample_id`, pattern = ":.+")) %>%
    dplyr::select(-`#sample_id`, -given_ancestry)
  somalier_all = dplyr::select(somalier_all, sample_id, predicted_ancestry, PC1, PC2, PC3, PC4, PC5)
  sample_table = left_join(sample_table, somalier_all)
  return(sample_table)
}



#' @title Collate ASHM Results.
#'
#' @description Determine the hypermutation status of a few genes.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr tidyr tibble
#'
#' @noRd
#'
#' @examples
#' sample_table = collate_ashm_results(sample_table = sample_table)
#'
#' @export
collate_ashm_results = function(sample_table,
                                seq_type_filter = "genome"){

  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter = seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  #just annotate BCL2, MYC and CCND1 hypermutation
  regions_df = data.frame(name = c("CCND1","BCL2","MYC"),
  region = c("chr11:69455000-69459900", "chr18:60983000-60989000", "chr8:128747615-128751834"))
  region_mafs = lapply(regions_df$region, function(x){get_ssm_by_region(region = x, streamlined = FALSE,seq_type=seq_type_filter)})
  tibbled_data = tibble(region_mafs, region_name = regions_df$name)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start, sample_id, region_name)

  tallied = unlisted_df %>%
    group_by(sample_id, region_name) %>%
    tally() %>%
    pivot_wider(values_from = n, names_from = region_name, values_fill = 0, names_prefix = "ashm_")

  sample_table = left_join(sample_table, tallied, by = "sample_id")
}



#' @title Collate CSR Results
#'
#' @description Collate a few CSR annotations, including MiXCR.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr
#'
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_csr_results(gambl_results_derived)
#'
#' @export
collate_csr_results = function(sample_table,
                               seq_type_filter = "genome"){

   if(seq_type_filter=="capture"){
     return(sample_table) #result doesn't exist or make sense for this seq_type
   }
   csr = suppressMessages(read_tsv("/projects/rmorin/projects/gambl-repos/gambl-nthomas/results/icgc_dart/mixcr_current/level_3/mixcr_genome_CSR_results.tsv"))
   sm_join = inner_join(sample_table, csr, by = c("sample_id" = "sample"))
   pt_join = inner_join(sample_table, csr, by = c("patient_id" = "sample"))
   complete_join = bind_rows(pt_join, sm_join) %>%
     bind_rows(dplyr::filter(sample_table, !patient_id %in% c(pt_join$patient_id, sm_join$patient_id))) %>%
     unique()

  return(complete_join)
}



#' @title Collate Curated SV Results.
#'
#' @description Collate all SV calls from the genome data and summarize for main oncogenes of interest per sample.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return The sample table with additional columns.
#'
#' @import readr dplyr
#'
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_curated_sv_results(gambl_results_derived)
#'
#' @export
collate_curated_sv_results = function(sample_table,
                                      seq_type_filter = "genome"){

  path_to_files = check_config_value(config::get("derived_and_curated"))
  project_base = check_config_value(config::get("project_base"))
  manual_files = dir(paste0(project_base, path_to_files), pattern = ".tsv")
  for(f in manual_files){
    full = paste0(project_base, path_to_files, f)
    this_data = suppressMessages(read_tsv(full, comment = "#"))
    #TO DO: fix this so it will join on biopsy_id or sample_id depending on which one is present, Done?
    sample_table = left_join(sample_table, this_data)
  }
  return(sample_table)
}



#' @title Collate Derived Results.
#'
#' @description Extract derived results stored in the database (these are usually slower to derive on the fly).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations. Default is FALSE, TRUE is not yet implemented.
#'
#' @return Data frame with one row per sample. Contains the contents of the derived_data table in the database.
#'
#' @import dplyr DBI RMariaDB
#'
#' @noRd
#'
#' @examples
#' gambl_results_derived = collate_derived_results(samples_df)
#'
#' @export
collate_derived_results = function(sample_table,
                                   seq_type_filter = "genome",
                                   from_flatfile = FALSE){

  if(from_flatfile){
    message("not implemented YET")
  }else{
    database_name = check_config_value(config::get("database_name"))
    con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = database_name)
    derived_tbl = dplyr::tbl(con, "derived_data") %>%
      as.data.frame()
  }
  derived_tbl = derived_tbl %>%
    dplyr::select(where( ~!all(is.na(.x)))) %>%
    as.data.frame() #drop the columns that are completely empty

  print(derived_tbl)
  sample_table = dplyr::left_join(sample_table, derived_tbl)
  print(sample_table)
  return(sample_table)
}



#' @title Collate Extra Metadata.
#'
#' @description Gather additional metadata information and expand the incoming sample table (or metadata).
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param file_path Path to extra metadata.
#'
#' @return A table.
#'
#' @import readr dplyr
#'
#' @noRd
#'
#' @examples
#' table = collate_extra_metadata(sample_table = "my_sample_table.txt")
#'
#' @export
collate_extra_metadata = function(sample_table,
                                  file_path){

  file_path = "/projects/rmorin/projects/gambl-repos/gambl-mutect2-lhilton/experiments/2021-04-21-Trios-MiXCR/trios_relapse_timing.tsv"
  extra_df = suppressMessages(read_tsv(file_path))
  sample_table = left_join(sample_table, extra_df, by = c("sample_id" = "biopsy_id"))
}



#' @title Collate NFKBIZ Results.
#'
#' @description Determine which cases have NFKBIZ UTR mutations.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#'
#' @return Samples table.
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' sample_table = collate_nfkbiz_results(sample_table = sample_table)
#'
#' @export
collate_nfkbiz_results = function(sample_table,
                                  seq_type_filter = "genome"){

  #TO DO: Update to work with hg38 projection
  if(missing(sample_table)){
    sample_table = get_gambl_metadata(seq_type_filter=seq_type_filter) %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  this_region = "chr3:101578214-101578365"
  nfkbiz_ssm = get_ssm_by_region(region = this_region,seq_type = seq_type_filter) %>%
    pull(Tumor_Sample_Barcode) %>%
    unique
  if(seq_type_filter=="genome"){
    nfkbiz_sv = get_manta_sv(region = this_region) %>%
      pull(tumour_sample_id) %>%
      unique
    nfkbiz = unique(c(nfkbiz_ssm, nfkbiz_sv))
  }
  else{
    nfkbiz = unique(nfkbiz_ssm)
  }


  sample_table$NFKBIZ_UTR = "NEG"
  sample_table[sample_table$sample_id %in% nfkbiz, "NFKBIZ_UTR"] = "POS"
  return(sample_table)
}



#' @title Collate PGA results for samples with CN data.
#'
#' @description Expand a metadata table horizontally with PGA metrics.
#'
#' @details Helper function called by `collate_results`, not meant for out-of-package usage.
#'
#' @param these_samples_metadata The metadata to be expanded with sample_id column.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @noRd
#'
#' @return data frame
#' @import dplyr
#'
#' @examples
#' # For genomes
#' meta <- get_gambl_metadata()
#' pga_metrics <- collate_pga(these_samples_metadata = meta)
#' # For exomes
#' meta_capture <- get_gambl_metadata(seq_type_filter = "capture")
#' pga_metrics_capture <- collate_pga(these_samples_metadata = meta_capture)
#'
#' @export
collate_pga <- function(
    these_samples_metadata,
    this_seq_type = "genome"
) {

    message(
        "Collating the PGA results ..."
    )
    # Currently only works for genomes
    if(! this_seq_type %in% c("genome", "capture")) {
        stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
    }

    # Default to all samples if sample table is missing
    if (missing(these_samples_metadata)) {
        message("No sample table was provided. Defaulting to all metadata ...")
        these_samples_metadata <- get_gambl_metadata(
            seq_type_filter = this_seq_type
        )
    }

    # Get the CN segments
    multi_sample_seg <- get_sample_cn_segments(
        sample_list = these_samples_metadata$sample_id,
        multiple_samples = TRUE,
        this_seq_type = this_seq_type
    ) %>%
    dplyr::rename("sample" = "ID")

    these_samples_pga <- calculate_pga(
        this_seg = multi_sample_seg
    )

    these_samples_metadata <- left_join(
        these_samples_metadata,
        these_samples_pga
    )

    return(these_samples_metadata)

}



#' @title Collate Quality Control Results.
#'
#' @description Expand a metadata table horizontally with quality control metrics.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table df with sample ids in the first column.
#' @param seq_type_filter default is genome, capture is also available for unix_group icgc_dart.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr readr
#'
#' @noRd
#'
#' @examples
#' qc_metrics = collate_qc_results(sample_table = sample_table)
#'
#' @export
collate_qc_results = function(sample_table,
                              seq_type_filter = "genome"){

  if(! seq_type_filter %in% c("genome", "capture")){
    stop("Please provide a valid seq_type (\"genome\" or \"capture\").")
  }

  #get paths
  base = check_config_value(config::get("project_base"))
  qc_template = check_config_value(config::get("qc_met"))

  #icgc_dart
  unix_group = "icgc_dart"
  icgc_qc_path = glue::glue(qc_template)
  icgc_qc_path_full = paste0(base, icgc_qc_path)

  #gambl
  unix_group = "gambl"
  gambl_qc_path = glue::glue(qc_template)
  gambl_qc_path_full = paste0(base, gambl_qc_path)

  #read in icgc qc data, rename sample id column and filter on samples in sample id in sample_table
  icgc_qc = suppressMessages(read_tsv(icgc_qc_path_full)) %>%
      dplyr::rename(sample_id = UID) %>%
      dplyr::select(-SeqType)

  #read in gambl qc data (if seq_type_filter set to "genome"), rename sample id column and filter on samples in sample id in sample_table
  if(seq_type_filter == "genome"){
    gambl_qc = suppressMessages(read_tsv(gambl_qc_path_full)) %>%
      dplyr::rename(sample_id = UID) %>%
      dplyr::select(-SeqType)

    #join gambl and icgc QC data
    full_qc = rbind(gambl_qc, icgc_qc)
    sample_table = left_join(sample_table, full_qc)

    #print n samples with QC metrics
    qc_samples = length(unique(full_qc$sample_id))
    message(paste("QC metrics for", qc_samples, "samples retrieved."))

  }else{
    message("Currently, seq_type_filter = \"capture\" is only available for unix_group \"icgc_dart\". Only QC metrics for icgc_dart will be returned.")
    #TO DO: Remove this once capture metrics have been generated for gambl samples.
    sample_table = left_join(sample_table, icgc_qc)

    #print n samples with QC metrics
    qc_samples = length(unique(icgc_qc$sample_id))
    message(paste("QC metrics for", qc_samples, "samples retrieved."))
  }
  return(sample_table)
}



#' @title Collate SBS Results.
#'
#' @description Bring in the results from mutational signature analysis.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param file_path Optional path to SBS file.
#' @param scale_vals Parameter not used?
#' @param sbs_manipulation Optional variable for transforming sbs values (e.g log, scale).
#'
#' @return A data frame with new columns added.
#'
#' @import dplyr tibble
#'
#' @noRd
#'
#' @examples
#' collated = collate_sbs_results(sample_table = sample_table,
#'                                sbs_manipulation = sbs_manipulation)
#
#' @export
collate_sbs_results = function(sample_table,
                               seq_type_filter = "genome",
                               file_path,
                               scale_vals = FALSE,
                               sbs_manipulation = ""){
  if(seq_type_filter!="genome"){
    message("skipping sbs for seq_type")
    return(sample_table)
  }
  if(missing(file_path)){
    base = check_config_value(config::get("project_base"))

    file_path = paste0(base,"icgc_dart/sigprofiler-1.0/02-extract/genome--hg38/BL_HGBL_DLBCL_FL_COMFL_CLL_MCL_B-ALL_PBL_DLBCL-BL-like_UNSPECIFIED_SCBC_MM_all/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt")
  }
  message(paste("loading",file_path))
  signatures = read.csv(file_path, sep = "\t", header = 1, row.names = 1)
  rs = rowSums(signatures)
  cn = colnames(signatures)
  new_sig = signatures
  if(sbs_manipulation == "scale"){
    #for(col in cn){
    #  scaled_vals = signatures[,col] / rs
    #  new_sig[,col]=scaled_vals
    #}
    #sbs1 = signatures[,"SBS1"] / rs
    #sbs5 = signatures[,"SBS5"] / rs
    #sbs9 = signatures[,"SBS9"] / rs
    #sbs8 = signatures[,"SBS8"] / rs
    #sbs12 = signatures[,"SBS12"] / rs
    #sbs17b = signatures[,"SBS17b"]/rs
    #sbs18 = signatures[,"SBS18"]/rs
    #sbs84 = signatures[,"SBS84"]/rs
    #sbs85 = signatures[,"SBS85"]/rs
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = apply(signatures, 2, function(x){x/rs}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
  }
  else if(sbs_manipulation == "log"){
    #sbs1 = log(signatures[,"SBS1"]+1)
    #sbs5 = log(signatures[,"SBS5"]+1)
    #sbs9 = log(signatures[,"SBS9"]+1)
    #sbs8 = log(signatures[,"SBS8"]+1)
    #sbs12 = log(signatures[,"SBS12"]+1)
    #sbs17b = log(signatures[,"SBS17b"]+1)
    #sbs18 = log(signatures[,"SBS18"]+1)
    #sbs84 = log(signatures[,"SBS84"]+1)
    #sbs85 = log(signatures[,"SBS85"]+1)
    sbs = apply(signatures, 2, function(x){log(x + 1)}) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
  }else{
    #sbs1 = signatures[,"SBS1"]
    #sbs5 = signatures[,"SBS5"]
    #sbs9 = signatures[,"SBS9"]
    #sbs8 = signatures[,"SBS8"]
    #sbs12 = signatures[,"SBS12"]
    #sbs17b = signatures[,"SBS17b"]
    #sbs18 = signatures[,"SBS18"]
    #sbs84 = signatures[,"SBS84"]
    #sbs85 = signatures[,"SBS85"]
    #sbs = data.frame(sample_id = rownames(signatures),sbs1=sbs1,sbs5=sbs5,sbs9=sbs9,sbs8=sbs8,sbs12=sbs12,sbs17b=sbs17b,sbs18=sbs18,sbs84=sbs84,sbs85=sbs85)
    sbs = signatures %>%
      rownames_to_column("sample_id")
  }
  sample_table = left_join(sample_table, sbs)
  return(sample_table)
}



#' @title Collate SSM Results.
#'
#' @description Compute summary statistics based on SSM calls.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param projection Specifies the projection, default is "grch37".
#' @param from_flatfile Optional argument whether to use database or flat file to retrieve mutations, default is TRUE.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return The sample table with additional columns.
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' ssm_results = colalte_ssm_results(sample_table = samples,
#'                                   include_silent = TRUE)
#'
#' @export
collate_ssm_results = function(sample_table,
                               seq_type_filter = "genome",
                               projection = "grch37",
                               from_flatfile = TRUE,
                               include_silent = FALSE){

  if(!include_silent){
    coding_class = coding_class[coding_class != "Silent"]
  }
  seq_type = seq_type_filter
  #iterate over every sample and compute some summary stats from its MAF
  if(from_flatfile){
    base_path = check_config_value(config::get("project_base"))
    #test if we have permissions for the full gambl + icgc merge
    maf_partial_path = check_config_value(config::get("results_flatfiles")$ssm$template$merged$deblacklisted)

    maf_path = paste0(base_path, maf_partial_path)
    maf_path = glue::glue(maf_path)
    message(paste("Checking permissions on:",maf_path))
    maf_permissions = file.access(maf_path, 4)
    if(maf_permissions == -1){
      message("fail. You do not have permissions to access all the results. Use the cached results instead.")
      return(sample_table)
    }
    print(paste("loading",maf_path))
    muts = vroom::vroom(maf_path) %>%
      dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification,t_alt_count,t_ref_count)
    mutated_samples = length(unique(muts$Tumor_Sample_Barcode))
    message(paste("mutations from", mutated_samples, "samples"))
  }
  #get tally of total per sample
  muts = muts %>%
    dplyr::rename("sample_id" = "Tumor_Sample_Barcode")

  muts = mutate(muts, vaf = t_alt_count/(t_alt_count + t_ref_count))
  muts_count = dplyr::select(muts, sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("total_ssm" = "n")

  sample_table = left_join(sample_table, muts_count)
  muts_mean = muts %>%
    dplyr::select(sample_id, vaf) %>%
    group_by(sample_id) %>%
    summarize(mean_vaf = mean(vaf))

  coding_mut = dplyr::filter(muts, Variant_Classification %in% coding_class)
  coding_mut_count = coding_mut %>%
    dplyr::select(sample_id) %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("coding_ssm" = "n")

  sample_table = left_join(sample_table, muts_mean)
  sample_table = left_join(sample_table, coding_mut_count)
  #check for coding SSMs in lymphoma genes
  coding_nhl = coding_mut %>%
    dplyr::filter(Hugo_Symbol %in% lymphoma_genes$Gene)

  coding_nhl_count = coding_nhl %>%
    group_by(sample_id) %>%
    tally() %>%
    dplyr::rename("driver_ssm" = "n")

  return(sample_table)
}



#' @title Collate SV Results.
#'
#' @description Determine and summarize which cases have specific oncogene SVs.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::collate_results], not meant for out-of-package usage.
#'
#' @param sample_table A data frame with sample_id as the first column.
#' @param tool Name of tool (optional, default is manta).
#' @param seq_type_filter Filtering criteria, default is genomes.
#' @param oncogenes Which oncogenes to collate SVs from.
#'
#' @return Data frame with additional columns ({tool}_{oncogene} and {tool}_{oncogene}_{partner}).
#'
#' @import dplyr
#'
#' @noRd
#'
#' @examples
#' results = collate_samples_sv_results(sample_table = samples,
#'                                      tool = "manta",
#'                                      oncogenes = c("MYC", "BCL2"))
#'
#' @export
collate_sv_results = function(sample_table,
                              tool = "manta",
                              seq_type_filter = "genome",
                              oncogenes = c("MYC", "BCL2", "BCL6", "CCND1", "IRF4")){
  if(seq_type_filter!="genome"){
    message("skipping sv for this seq_type")
    return(sample_table)
  }
  if(tool == "manta"){
    all_svs = get_manta_sv()
  }
  annotated_svs = annotate_sv(all_svs) %>%
  dplyr::filter(!is.na(partner))

  if(missing(sample_table)){
    sample_table = get_gambl_metadata() %>%
      dplyr::select(sample_id, patient_id, biopsy_id)
  }
  multiout = function(df,
                       annotated,
                       tool,
                       oncogene_name){

    some_fusions = dplyr::filter(annotated, gene == all_of(oncogene_name)) %>%
      group_by(tumour_sample_id) %>%
      arrange(partner) %>%
      dplyr::filter(row_number() == 1)

    df = mutate(df, "{tool}_{oncogene_name}_sv" := case_when(sample_id %in% some_fusions$tumour_sample_id ~ "POS", TRUE ~ "NEG"))
    some_fusions = some_fusions %>%
      dplyr::select(tumour_sample_id, partner) %>%
      mutate("{tool}_{oncogene_name}_partner" := partner) %>%
      dplyr::select(-partner)

    df = left_join(df, some_fusions, by = c("sample_id" = "tumour_sample_id"))
    return(df)
  }
  out_table = sample_table
  for(oncogene in oncogenes){
    out_table = multiout(out_table, annotated_svs, "manta", oncogene)
  }
  return(out_table)
}



#' @export
compare_coding_mutation_pattern = function(maf_df1,maf_df2,gene){
  if(missing(maf_df1) | missing(maf_df2)){
    stop("must provide two data frames containing mutations you would like to compare")
  }
  if(missing(gene)){
    stop("Must provide the Hugo_Symbol of a single gene that is present in both maf files")
  }
  missense_positions1 = dplyr::filter(maf_df1,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% str_remove_all("p.\\w") %>% str_extract("\\d+") %>% as.numeric()
  missense_positions2 = dplyr::filter(maf_df2,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% str_remove_all("p.\\w") %>% str_extract("\\d+") %>% as.numeric()
 if(length(missense_positions1)==0 | length(missense_positions2)==0 ){
   message(paste("no mutations for",gene,"in one or both data sets"))
   return(list(kl=15))
 }
  #generate range of amino acids based on what we can infer from the MAF (not ideal)
  max_pos = max(c(missense_positions1,missense_positions2))
  full_df = data.frame(position=c(1:max_pos))
  df1 = data.frame(position=missense_positions1) %>% group_by(position) %>% tally() %>% rename("group1"="n")
  df2 = data.frame(position=missense_positions2) %>% group_by(position) %>% tally() %>% rename("group2"="n")
  full_df = left_join(full_df,df1) %>% mutate(group1=ifelse(is.na(group1),0,group1))
  full_df = left_join(full_df,df2) %>% mutate(group2=ifelse(is.na(group2),0,group2))
  # convert to the format needed by KL
  all_counts = dplyr::select(full_df,-position) %>% t()
  all_counts[1,]=all_counts[1,]/sum(all_counts[1,])
  all_counts[2,]=all_counts[2,]/sum(all_counts[2,])
  kl_out = KL(all_counts)
  return(list(df=full_df,kl=unname(kl_out)))
}



#' @title Compare Mutation Flavour.
#'
#' @description Get a MAF that is just the variants unique to one of two flavours of variant calls available.
#'
#' @details Subset a MAF to only have variants that are unique to one flavour (specified with `flavour1`).
#' This function is currently not exported, since there is only one flavour available at the moment (see docs for [GAMBLR::get_ssm_by_sample]).
#'
#' @param these_sample_ids A vector of sample IDs to be included.
#' @param flavour1 First flavour of variant calls, to be returned as unique if not present in flavour2. Default is "clustered".
#' @param flavour2 Second flavour of variant calls.
#'
#' @return a list with MAFs that are only present in flavour1.
#'
#' @noRd
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @export
compare_mutation_flavour = function(these_sample_ids,
                                    flavour1 = "clustered",
                                    flavour2 = ""){

  these_dfs = list()
  for(this_sample_id in these_sample_ids){
    message(this_sample_id)
    maf1 = get_ssm_by_sample(this_sample_id, flavour = flavour1, this_seq_type = seq_type)
    maf2  = get_ssm_by_sample(this_sample_id, flavour = flavour2)
    maf1_only = intersect_maf(maf1, maf2)
    these_dfs[[this_sample_id]] = maf1_only
  }
  this_maf = rbindlist(these_dfs, use.names = TRUE)
  return(this_maf)
}



#' @title Refresh Metadata Tables
#'
#' @description Refresh the contents of a metadata table.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package usage.
#'
#' @return Table.
#'
#' @import RMariaDB DBI dplyr
#'
#' @noRd
#'
#' @examples
#' ref_meta = referesh_metadata_tables()
#'
#' @export
referesh_metadata_tables = function(){

  con = dbConnect(RMariaDB::MariaDB(), dbname = database_name)
  all_metadata_info = sanity_check_metadata()
  tables = pull(all_metadata_info, table)
  files = pull(all_metadata_info, file)
  for(i in c(1:length(files))){
    refresh_full_table(tables[i], con, files[i])
  }
}



#' @title Refresh Full Table
#'
#' @description Refresh the contents of a database table.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::referesh_metadata_tables], not meant for out-of-package usage.
#'
#' @param table_name Name of table to refresh.
#' @param connection Database connection object.
#' @param file_path Path to the table contents to populate.
#'
#' @return A table.
#'
#' @import DBI RMariaDB readr
#'
#' @noRd
#'
#' @examples
#' refresh_full_table(table_x, con,file_x)
#'
#' @export
refresh_full_table = function(table_name,
                              connection,
                              file_path){

  table_data = suppressMessages(read_tsv(file_path))
  dbWriteTable(con, table_name, table_data, overwrite = TRUE)
  print(paste("POPULATING table:", table_name, "USING path:", file_path))
}



#' @title Region To Chunks.
#'
#' @description Parse a region string into; chromosome, start and end.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::calc_mutation_frequency_sliding_windows], not meant for out-of-package usage.
#'
#' @param region A region string e.g. "chrX:12345-678910".
#'
#' @return A named list.
#'
#' @noRd
#'
#' @examples
#' chr_start_end = region_to_chunks("chr1:1111-2222")
#'
#' @export
region_to_chunks = function(region){

  region = unname(region)
  region = gsub(",", "", region)
  #format is chr6:37060224-37151701
  split_chunks = unlist(strsplit(region, ":"))
  chromosome = split_chunks[1]
  startend = unlist(strsplit(split_chunks[2], "-"))
  qstart = startend[1]
  qend = startend[2]
  return(list(chromosome = chromosome, start = qstart, end = qend))
}



#' @title Sanity Check Metadata.
#'
#' @description Function that performs sanity checks on metadata.
#'
#' @details Helper function for sanity checking GAMBL metadata.
#'
#' @return A table.
#'
#' @import tibble readr dplyr
#'
#' @noRd
#'
#' @examples
#' sane_meta_data = sanity_check_metadata()
#'
#' @export
sanity_check_metadata = function(){

  cfg = check_config_value(config::get("tables"))
  database_name = check_config_value(config::get("database_name"))
  metadata_tables = tibble(key = names(cfg), table = cfg) %>%
    unnest_auto("table")

  cfg = check_config_value(config::get("table_flatfiles"))
  metadata_files = tibble(key = names(cfg), file = cfg) %>%
    unnest_auto("file")

  all_metadata_info = left_join(metadata_tables, metadata_files)
  base_path = check_config_value(config::get("repo_base"))
  all_metadata_info = all_metadata_info %>%
    mutate(file = paste0(base_path, file))

  all_metadata_df = all_metadata_info %>%
    column_to_rownames(var = "key")
  #all samples with different seq_type and protocol must have a unique sample_id
  sample_df = suppressMessages(read_tsv(all_metadata_df["samples", "file"]))
  tumour_samples = sample_df %>%
    dplyr::select(patient_id, sample_id, biopsy_id, seq_type, protocol) %>%
    dplyr::filter(!is.na(biopsy_id))

  n_samp_bio = tumour_samples %>%
    count() %>%
    pull(n)

  #check for any multiplicity of sample_id
  n_samp = tumour_samples %>%
    dplyr::select(-biopsy_id) %>%
    unique() %>%
    count() %>%
    pull(n)

  #should be the same number as above
  if(!n_samp == n_samp_bio){
    print("ERROR! some biopsies appear to have the same sample_id/protocol combination")
  }
  return(all_metadata_info)
}



#' @export
socketWrite = function (sock, string) {
  print(string)
  write.socket(sock, string)
  response <- read.socket(sock)
  return(response)
}



#' @title Standardize Chromosome Prefix.
#'
#' @description Standardize the chr prefix in a vector of chromosome names based on projection.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package use.
#'
#' @param incoming_vector Input vector of any length with chromosome names.
#' @param projection Projection to which chr prefix should be standardized.
#'
#' @return A vector of chromosome names with prefix standardized to projection
#'
#' @noRd
#'
#' @examples
#' these_chrs = c(8, "13", "chr4", "chrY")
#'
#' standardize_chr_prefix(incoming_vector = these_chrs,
#'                        projection = "hg38")
#'
#' @export
standardize_chr_prefix = function(incoming_vector,
                                  projection){

  if (projection %in% c("grch37", "grch38")) {
    output_vector = gsub("chr", "", incoming_vector)
  } else {
    output_vector = gsub("chr", "", incoming_vector) # if there is a mix of prefixed and non-prefixed options
    output_vector = paste0("chr", output_vector)
  }
  return(output_vector)
}



#' @title Subset CN States.
#'
#' @description Get the available CN states in the incoming data frame.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR::fancy_multisample_ideo], for sub-setting copy number information based on segments available in cn data
#'
#' @param cn_segments DF with copy number segments, usually retrieved from get_sample_cn_segments.
#' @param include_2 Optional parameter for including or omit CN state == 2. Default is FALSE.
#' @param samplen Numeric value that annotates the sample order.
#'
#' @return Nothing.
#'
#' @noRd
#'
#' @examples
#' cn_states = get_sample_cn_segments(multiple_samples = TRUE,
#'                                    sample_list = c("00-15201_tumorA",
#'                                                    "HTMCP-01-06-00422-01A-01D"),
#'                                    streamlined = FALSE)
#'
#' subset_cnstates(cn_segments = cn_states,
#'                 samplen = 1)
#'
#' @export
subset_cnstates = function(cn_segments,
                           include_2 = FALSE,
                           samplen){

  #transform CN states > 6 = 6 (to reflect the current copy number palette in gamblr)
  cn_segments$CN[cn_segments$CN > 6] = 6

  #filter out CN == 2
  if(!include_2){
    cn_segments = subset(cn_segments, CN != 2)
  }

  #update CN annotations (if present in cn_segment data).
  cn_segments$CN = paste0("cn_", cn_segments$CN , "_sample", samplen)

  #convert to factor.
  cn_segments$CN = as.factor(cn_segments$CN)

  #split cn_segments on available factors and lists into the global environment.
  l = split(cn_segments, cn_segments$CN)
  list2env(l, envir = .GlobalEnv)
}



#' @export
test_glue = function(placeholder="INSERTED"){
  some_string = "this text has {placeholder}"
  print(some_string)
  ss=glue::glue(some_string)
  print(ss)
}
