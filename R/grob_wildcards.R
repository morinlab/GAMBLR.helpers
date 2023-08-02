#Helper functions not for export
#' @export
grob_wildcards = function(wildcarded_string){
  wildcards = unlist(stringr::str_extract_all(wildcarded_string,"\\{[^\\{]+\\}"))
  wildcards = stringr::str_remove_all(wildcards,"\\{") %>%  stringr::str_remove_all(.,"\\}")
  return(wildcards)
}
