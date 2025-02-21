#' @title Check Config Value.
#'
#' @description Check the existence of a specific config key.
#' The function will notify the user and end the program if no such key exists.
#'
#' @details INTERNAL FUNCTION for checking the existence of a config value,
#' not meant for out-of-package usage.
#'
#' @param config_key key from config, prefixed with config::get()
#'
#' @return A string with the path to a file specified in the config or nothing
#' (if config key is NULL).
#' 
#' @examples
#' check_config_value(config::get("resources")$blacklist$template)
#'
#' @export
check_config_value = function(config_key){
  if(is.null(config_key)){
    stop(paste0("ATTENTION! The above described key is missing from the config, make sure your config is up to date"))
  }else{
    return(config_key)
  }
}

#' @title Check Config and Value.
#'
#' @description Check the existence of a specific config key.
#' The function will notify the user and end the program if no such key exists.
#'
#' @details INTERNAL FUNCTION for checking the existence of a config value,
#' not meant for out-of-package usage.
#'
#' @param config_key A string of one or more nested keys to retrieve 
#' from the config in the format "$key1$key2$keyN"
#'
#' @return A string with the value from the config key or nothing
#' (if config key is NULL or missing).
#' 
#' @examples
#' check_config_value(config::get("resources")$blacklist$template)
#'
#' @export
check_config_and_value <- function(config_key) {
  # Check if config.yml is in working directory
  if (file.exists("config.yml")) {
    config_list <- config::get()
  } else {
    # If not, fall back to package-specific config
    print("config.yml not in working directory")
    calling_pkg <- get_calling_package()
    print(paste("function called from:", calling_pkg))
    if(is.null(calling_pkg)){
      calling_pkg <- "GAMBLR.results"
    }
    # Check for package-specific config file
    config_file <- system.file("extdata", "config.yml", package = calling_pkg)
    
    if (file.exists(config_file) && nzchar(config_file)) {
      message("Using package-specific config file: ", config_file)
      config_list <- config::get(file = config_file)
    } else {
      print(config_file)
      stop("No config.yml found in working directory or package. Check installation and file paths.")
    }
  }
  
  # Navigate the config list to find the requested key
  if (!is.null(config_key)) {
    
    # Split the key by "$" and remove empty parts
    key_parts <- unlist(strsplit(config_key, "\\$"))
    key_parts <- key_parts[key_parts != ""]  # Remove empty parts
    
    # Recursively traverse the list to find the value
    get_nested_value <- function(list_obj, keys) {
      if (length(keys) == 0) {
        return(list_obj)
      } else {
        first_key <- keys[1]
        remaining_keys <- keys[-1]
        
        if (is.list(list_obj) && first_key %in% names(list_obj)) {
          return(get_nested_value(list_obj[[first_key]], remaining_keys))
        } else {
          stop(paste0("Config key not found: ", config_key))
        }
      }
    }
    
    # Start the recursive search
    value <- get_nested_value(config_list, key_parts)
    return(value)
  } else {
    stop("No config key provided.")
  }
}





get_calling_package <- function() {
  # Get the call stack and parents
  calls <- sys.calls()
  parents <- sys.parents()
  
  # Loop through the call stack from the oldest to the newest call
  for (i in seq_along(parents)) {
    # Get the environment of the parent frame
    env <- sys.frame(i)
    
    # Get the environment one level up
    parent_env <- parent.env(env)
    
    # If it's a namespace environment, get the package name
    if (isNamespace(parent_env)) {
      pkg_name <- getNamespaceName(parent_env)
      
      # Skip the current package's namespace
      if (pkg_name != "GAMBLR.helpers") {  # Change this to your current package's name
        return(pkg_name)
      }
    }
  }
  
  # If no package is found, return NULL
  return(NULL)
}