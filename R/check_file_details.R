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
