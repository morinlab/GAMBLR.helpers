#' @title Glob wildcards.
#'
#' @description Helper function to extract wildcard values from a string.
#'
#' @param wildcarded_string String containing wildcards inside the {} notations.
#'
#' @return vector
#'
#' @export
grob_wildcards = function(wildcarded_string){
  wildcards = unlist(regmatches(wildcarded_string, gregexpr("\\{[^\\{]+\\}", wildcarded_string)))
  wildcards = gsub("\\{", "", wildcards) %>% gsub("\\}", "", .)
  return(wildcards)
}
