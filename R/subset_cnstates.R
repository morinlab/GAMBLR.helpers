#' @title Subset CN States.
#'
#' @description Get the available CN states in the incoming data frame.
#'
#' @details For sub-setting copy number information based on segments available in cn data
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
#' cn_states = get_sample_cn_segments(these_sample_ids = c("00-15201_tumorA",
#'                                                         "HTMCP-01-06-00422-01A-01D"),
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
