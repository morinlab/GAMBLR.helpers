#' @title Check Config Value.
#'
#' @description Check the existence of a specific config key.
#' The function will notify the user and end the program if no such key exists.
#'
#' @details INTERNAL FUNCTION for checking the existence of a config value, not meant for out-of-package usage.
#'
#' @param config_key key from config, prefixed with config::get()
#'
#' @return A string with the path to a file specified in the config or nothing (if config key is NULL).
#' 
#' @noRd
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
