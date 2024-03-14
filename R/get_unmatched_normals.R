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
#' @import dplyr tidyr tibble
#' @export
get_unmatched_normals = function(seq_type_filter){
    a <- check_config_value(config::get("unmatched_normal_ids"))
    df <- a %>%
        as.data.frame(
            check.names = FALSE
        ) %>%
        t %>%
        as.data.frame() %>%
        tibble::rownames_to_column("rownames") %>%
        tidyr::separate(
            rownames,
            into = c("unix_group", "seq_type", "genome_build"),
            sep = "\\."
        ) %>%
        dplyr::select(
            "normal_sample_id" = "V1",
            genome_build, seq_type, unix_group
        ) %>%
        dplyr::filter(seq_type == seq_type_filter)
    return(df)
}
