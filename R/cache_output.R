#' @title Cache output.
#'
#' @description TODO.
#'
#' @param result_df TODO.
#' @param function_name TODO.
#' @param clobber_mode TODO.
#' @param get_existing TODO.
#' @param function_params TODO.
#' @param additional_details TODO.
#'
#' @return full path to the file that was written
#'
#' @export
cache_output = function(result_df,
                        function_name,
                        clobber_mode = F,
                        get_existing = F,
                        function_params = list(region = "chr3:98300000-198022430", bin_size=2000, seq_type="genome"),
                        additional_details = list(foreground = "DLBCL_FL_BL", background = "CLL_MM_MCL")){

  cache_file_name = paste0(check_config_value(config::get("repo_base")),"cached_results/", function_name)
  for (param in names(function_params)[order(names(function_params))]){
    cache_file_name = paste0(cache_file_name,"--",param,"-",function_params[[param]])
  }
  for (detail in names(additional_details)){
    cache_file_name = paste0(cache_file_name,"--",detail,"-",additional_details[[detail]])
  }
  cache_file_name = paste0(cache_file_name,".tsv")
  if(file.exists(cache_file_name)){
    if(get_existing){
      result_df = suppressMessages(read_tsv(cache_file_name))
      return(result_df)
    }
    if(!clobber_mode){
      warning(paste("file",cache_file_name,"exists!"))
      stop("Will not overwrite unless you rerun this in clobber_mode = TRUE")
    }
  }else{
    if(get_existing){
      stop(paste("cannot find cached result for this parameter combination",cache_file_name))
    }
  }

  message(paste("creating/overwriting",cache_file_name))
  write_tsv(result_df,file=cache_file_name)
  return(cache_file_name)
}
