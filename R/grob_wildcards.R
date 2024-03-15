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
  wildcards = unlist(stringr::str_extract_all(wildcarded_string,"\\{[^\\{]+\\}"))
  wildcards = stringr::str_remove_all(wildcards,"\\{") %>%  stringr::str_remove_all(.,"\\}")
  return(wildcards)
}
