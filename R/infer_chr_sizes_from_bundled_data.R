#' @title Infer chromosome sizes using the bundled data
#'
#' @param projection Specify which genome build to use. One of "grch37" (default) or "hg38".
#'
#' @return A data frame with columns `chromosome` and `size.`
#' 
#' @import dplyr GAMBLR.data
#' @export
#'
#' @examples
#' infer_chr_sizes_from_bundled_data(projection = "hg38")
#' 
infer_chr_sizes_from_bundled_data <- function(projection = "grch37"){
  if(projection == "grch37"){
    chr_arms <- GAMBLR.data::chromosome_arms_grch37
  }else if(projection == "hg38"){
    chr_arms <- GAMBLR.data::chromosome_arms_hg38
  }else{
    stop("projection parameter must be \"grch37\" or \"hg38\".")
  }
  chr_end <- dplyr::filter(chr_arms, arm == "q") %>% 
    pull(end) 
  dplyr::filter(chr_arms, arm == "p") %>% 
    mutate(size = chr_end - start + 1) %>% 
    dplyr::select(chromosome, size)
}
