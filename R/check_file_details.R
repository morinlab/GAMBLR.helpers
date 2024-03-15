#' @title Check file details
#'
#' @description When relative path to a file is given, this function will
#' automatically generate it's full path using the project_base config value
#' and check for existence of the file at the full path. Can operate on vector
#' of several relative paths.
#'
#' @param relative_paths Vector of relative paths.
#' @return Boolean
#'
#' @export
check_file_details = function(relative_paths){
  not_found = c()
  base_path = check_config_value(config::get("project_base"))
  for(relative_path in relative_paths){
    full_path = file.path(base_path, relative_path)
    message(paste("Looking for:",full_path))
    if(file.exists(full_path)){
      message("OK")
    }else{
      #message("Uh oh. This file cannot be found!")
      print(full_path)
      not_found=c(not_found,full_path)
      #print("-=10101010101=-")
    }
  }
  return(not_found)
}
