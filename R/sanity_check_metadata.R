#' @title Sanity Check Metadata.
#'
#' @description Function that performs sanity checks on metadata.
#'
#' @details Helper function for sanity checking GAMBL metadata.
#'
#' @return A table.
#'
#' @import tibble readr dplyr tidyr
#'
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
