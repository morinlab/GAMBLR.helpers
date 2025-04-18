#' @title Get template wildcards.
#'
#' @description TODO.
#'
#' @param parent_key TODO.
#' @param template_key TODO.
#'
#' @return list
#'
#' @export
get_template_wildcards = function(parent_key,
                                  template_key){

  if(missing(template_key)){
    wildcard_string = config::get(parent_key)
  }else{
    wildcard_string = config::get(paste0(parent_key,"_wildcards"))[template_key]
  }
  wildcards = strsplit(wildcard_string,",")
  return(unlist(wildcards))
}
