#' @title Standardize Chromosome Prefix.
#'
#' @description Standardize the chr prefix in a vector of chromosome names based
#' on projection.
#'
#' @details INTERNAL FUNCTION, not meant for out-of-package use.
#'
#' @param incoming_vector Input vector of any length with chromosome names.
#' @param projection Projection to which chr prefix should be standardized.
#'
#' @return A vector of chromosome names with prefix standardized to projection
#'
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
    # if there is a mix of prefixed and non-prefixed options
    output_vector = gsub("chr", "", incoming_vector)
    output_vector = paste0("chr", output_vector)
  }
  return(output_vector)
}
