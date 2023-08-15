#' @title Region To Chunks.
#'
#' @description Parse a region string into; chromosome, start and end.
#'
#' @details INTERNAL FUNCTION called by [GAMBLR.utils::calc_mutation_frequency_bin_region], not meant for out-of-package usage.
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
