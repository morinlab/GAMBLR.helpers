# adding coding class to global environment
coding_class = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")



#' @title Populate Tool Results.
#'
#' @description Populate the database with the per-sample summarized results of various tools.
#'
#' @details this function is still in draft mode, export to NAMESPACE has been removed for now.
#'
#' @param tool_name Name of the tool to get the results for.
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' results = populate_tool_results(tool_name = "slims_3")
#' }
#' 
populate_tool_results = function(tool_name){

  #IMPORTANT TODO: This function should only ever work with samples that exist in the metadata
  # Perhaps it should report any excluded outputs in case they need to be deleted from the main output directories
  matched_analyses = unlist(check_config_value(config::get("analyses")$matched))
  print(matched_analyses)
  database_name = check_config_value(config::get("database_name"))
  genome_builds = unlist(strsplit(check_config_value(config::get("genome_builds")), ","))
  groups = unlist(strsplit(check_config_value(config::get("unix_groups")), ","))
  for(analysis_type in names(matched_analyses)){
    tool_name = matched_analyses[analysis_type]
    message(paste("populating results for", tool_name))
    populate_each_tool_result(tool = tool_name, genome_builds, groups)
  }
}



#' @title Process All Manta Bedpe.
#'
#' @description This function is in draft mode.
#'
#' @details This is a helper function that is not meant to be used routinely.
#'
#' @param file_df Paths to bedpe.
#' @param out_dir output directory.
#' @param group The unix group.
#' @param genome_build Genome build.
#' @param projection_build The genome we want all results to be relative to (lifted if necessary).
#'
#' @import dplyr readr
#' 
#' @noRd
#'
process_all_manta_bedpe = function(file_df,
                                   out_dir,
                                   group,
                                   genome_build,
                                   projection_build = "grch37"){

  to_merge = list()
  if(missing(out_dir)){
    project_base = check_config_value(config::get("project_base"))
    base_out_dir = check_config_value(config::get("results_staging")$manta)
    out_dir = paste0(project_base, group, "/", base_out_dir)
  }

  process_manta = function(bedpe_file, liftover_to_hg19 = FALSE, liftover_to_hg38 = FALSE, only_return_missing = FALSE, projection = "grch37"){
    cnames = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SOMATIC_SCORE", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "VAF_tumour", "VAF_normal", "DP_tumour", "DP_normal", "tumour_sample_id", "normal_sample_id", "pair_status")
    svbed = suppressMessages(read_tsv(bedpe_file, comment = "##", col_types = "cddcddccccccccccccccccc"))
    this_patient = colnames(svbed)[23]
    this_normal = colnames(svbed)[22]

    if(grepl("--unmatched", bedpe_file)){
      pair_status = "unmatched"
      svbed$pair_status = "unmatched"
    }else{
      svbed$pair_status = "matched"
      pair_status = "matched"
    }
    is_lifted = "native"
    if(liftover_to_hg19 || liftover_to_hg38){
      is_lifted = "lifted"
    }

    if(genome_build == projection | (genome_build == "hs37d5" & projection == "grch37")){
      is_lifted = "native"
    }
    out_file = paste0(out_dir, "/", this_patient, "--", this_normal, "--", pair_status, "--", is_lifted, "--genome--", genome_build, "--", projection_build, "_sv.tsv")
    message("working on OVER HERE:", bedpe_file)
    print(paste("output:", out_file))
    if(file.exists(out_file)){
      if(!only_return_missing){
        print(paste("LOADING", out_file))
        svbed = suppressMessages(read_tsv(out_file, col_types = "ccccccccccccnnnnccc", col_names = cnames))
        return(svbed)
      }
      else{
        svbed = dplyr::filter(svbed, is.na(tumour_sample_id))
        return(svbed)
      }
    }
    if(liftover_to_hg19){
      svbed = liftover_bedpe(bedpe_df = svbed)
    }else if(liftover_to_hg38){
      svbed = liftover_bedpe(bedpe_df = svbed, target_build = "hg38")
    }

    infos = pull(svbed, this_patient)
    infos_n = pull(svbed, this_normal)
    colnames(svbed)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")

    svbed$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
    svbed$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
    svbed$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
    svbed$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
    svbed$SOMATIC_SCORE = sapply(svbed$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

    svbed$tumour_sample_id = this_patient
    svbed$normal_sample_id = this_normal
    message(paste("checking status:", bedpe_file))

    svbed$NAME = "."
    svbed = svbed %>%
      dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, NAME, SOMATIC_SCORE, STRAND_A, STRAND_B, TYPE, FILTER, VAF_tumour, VAF_normal, DP_tumour, DP_normal, tumour_sample_id, normal_sample_id, pair_status)

    #remove chr prefix from both chromosome names
    svbed = svbed %>%
      mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
      mutate(CHROM_B = gsub("chr", "", CHROM_B))
    #print(paste("writing output to",out_file))
    #run liftover after formatting?

    write_tsv(svbed, out_file, col_names = FALSE)
    return(svbed)
  }

  #separately run the hg38 and other builds, separately run per unix_group
  if(projection_build == "grch37"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df, genome_build == "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_lifted = hg38_files %>%
        purrr::map(process_manta, liftover_to_hg19 = TRUE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df, genome_build != "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_not_lifted = not_hg38_files %>%
        purrr::map(process_manta, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }
  }else if(projection_build == "hg38"){
    if(genome_build == "hg38"){
      hg38_files = dplyr::filter(file_df, genome_build == "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_not_lifted = hg38_files %>%
        purrr::map(process_manta, liftover_to_hg38 =FALSE, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }else{
      not_hg38_files = dplyr::filter(file_df, genome_build != "hg38" & unix_group == group) %>%
        pull(file_path)

      bed_data_lifted = not_hg38_files %>%
        purrr::map(process_manta, liftover_to_hg38 = TRUE, liftover_to_hg19 = FALSE, projection = projection_build) %>%
        purrr::reduce(rbind)
    }
  }
}



#' @title Read Merge Manta With Liftover.
#'
#' @description Takes a path to bedpe and runs liftover ([GAMBLR::liftover_bedpe]) based on the original genome build of the bedpe.
#'
#' @details This is a helper function that is not meant to be used routinely.
#'
#' @param bedpe_paths path to bedpe
#' @param pattern pattern
#' @param out_dir output directory
#'
#' @import dplyr readr
#' 
#' @noRd
#'
#' @examples
#' \dontrun{
#' manta_bedpe = read_merge_manta_with_liftover(bedpe_paths = "some_path.bedpe",
#'                                              out_dir = "../")
#' }
#' 
read_merge_manta_with_liftover = function(bedpe_paths = c(),
                                          pattern = "--matched",
                                          out_dir){

  to_merge = list()
  print(head(bedpe_paths))
  for(thispath in bedpe_paths){
    sample_paths = dir(thispath, pattern = pattern) #DEBUGGING
    print(sample_paths)
    #sample_paths = head(sample_paths,15) #for debugging
    for(sp in sample_paths){
      full_path = paste0(thispath, sp, "/somaticSV.bedpe")
      print(paste("working on HERE:", full_path))
      if(grepl("hg38", full_path)){
        print("using liftOver")
        svbed = liftover_bedpe(full_path) #load and convert to grch37 coordinates
      }else{
        svbed = suppressMessages(read_tsv(full_path, comment = "##", col_types = "cddcddccccccccccccccccc"))
      }

      this_patient = colnames(svbed)[23]
      this_normal = colnames(svbed)[22]
      out_file = paste0(out_dir, "/", this_patient, "--", this_normal, "--hg38Togrch37_sv.tsv")
      print(paste("writing output to", out_file))
      infos = pull(svbed, this_patient)
      infos_n = pull(svbed, this_normal)
      colnames(svbed)[c(1:6)] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B")
      #all_vafs = get.sv.vaf(infos)
      #svbed$VAF = as.numeric(all_vafs)
      svbed$VAF_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
      svbed$DP_tumour = sapply(infos, function(x){as.numeric(tail(unlist(strsplit(x, ":")),2)[1])})
      svbed$VAF_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 1))})
      svbed$DP_normal = sapply(infos_n, function(x){as.numeric(tail(unlist(strsplit(x, ":")), 2)[1])})
      svbed$SOMATIC_SCORE = sapply(svbed$INFO_A, function(x){as.numeric(tail(unlist(strsplit(x, "=")), 1))})

      #filter on PASS, score, VAF
      #svbed_filt = svbed %>%
      #  filter(SCORE > minScore & FILTER == "PASS") %>%
      #  dplyr::select(c(chrom1, start1, end1, chrom2, start2, end2))

      svbed$tumour_sample_id = this_patient
      svbed$normal_sample_id = this_normal
      if(grepl("--unmatched", sp)){
        svbed$pair_status = "unmatched"
      }else{
        svbed$pair_status = "matched"
      }
      print(head(svbed))
      svbed$NAME = "."

      svbed = svbed %>%
        dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, NAME, SOMATIC_SCORE, STRAND_A, STRAND_B, TYPE, FILTER, VAF_tumour, VAF_normal, DP_tumour, DP_normal, tumour_sample_id, normal_sample_id, pair_status)
      #remove chr prefix from both chromosome names
      svbed = svbed %>%
        mutate(CHROM_A = gsub("chr", "", CHROM_A)) %>%
        mutate(CHROM_B = gsub("chr", "",CHROM_B))

      write_tsv(svbed, out_file, col_names = FALSE)
      #to_merge[[this_patient]] = svbed
    }
  }
}
