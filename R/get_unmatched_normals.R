#' @title Get unmatched normals.
#'
#' @description Helper function to get the unmatched normals from the main
#' config.
#'
#' @param seq_type_filter Seq type key from config for which to return the
#'      unmatched normals for.
#'
#' @return data frame
#'
#' @export
get_unmatched_normals = function(seq_type_filter){
  a = check_config_value(config::get("unmatched_normal_ids"))
  df = melt(a,value.name="normal_sample_id") %>%
    rename(c("genome_build"="L3","seq_type"="L2","unix_group"="L1")) %>%
    dplyr::filter(seq_type == seq_type_filter)
  return(df)
}
