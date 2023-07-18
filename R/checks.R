#Global variable specifying what metadata columns are absolutely required
#' @export
required_cols = c("sample_id","patient_id","pathology","seq_type","genome_build","pairing_status","Tumor_Sample_Barcode")



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



#' @export
check_expected_outputs = function(tool_name="battenberg",seq_type_filter="genome"){
  projection = get_template_wildcards("projections")
  #drop irrelevant rows of the metadata based on the scope of the tool etc
  if(tool_name=="battenberg"){
    template_path = check_config_value(config::get("results_flatfiles")$cnv$battenberg)
    extra_wildcards = check_config_value(config::get("results_flatfiles")$cnv$battenberg_wildcards)
    #in the current setup, this drops unmatched samples (could be hard-coded but using the config is more transparent)
    relevant_metadata = dplyr::filter(all_meta,base::get(names(extra_wildcards)[1]) == unname(extra_wildcards[1]))
    runs = dplyr::filter(relevant_metadata,seq_type == seq_type_filter) %>%
      mutate(tumour_sample_id = sample_id)
  }else if(tool_name=="slms_3"){
    runs = get_runs_table(seq_type_filter = seq_type_filter) %>% mutate(pair_status = pairing_status)
    vcf_base = get_template_wildcards("vcf_base_name")
    runs = mutate(runs,vcf_base_name = vcf_base)
    #runs = mutate(runs,target_builds = projection)
    runs = expand_grid(runs,target_builds=projection)
    template_path = check_config_value(config::get("results_flatfiles")$ssm$template$clustered$deblacklisted)
  }
  w = grob_wildcards(template_path)

  #use the unix group and seq_type from the actual metadata and all available projections
  seq_type = seq_type_filter

  runs_files = mutate(runs,outfile=glue::glue(template_path))

}



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



#' @export
check_times = function(relative_paths,
                       archive_mode = FALSE,
                       force_backup = FALSE){

  local_base = check_config_value(config::get("project_base"))
  remote_base = check_config_value(config::get("project_base",config="default"))
  if(archive_mode){
    if(local_base == remote_base){
      message("checking against local archive")
    }else{
      message("Currently, this mode must be run on the GSC (not remotely)")
      return(NULL)
    }
    local_base = check_config_value(config::get("archive"))

  }
  for(rel_f in relative_paths){
    local_f = paste0(local_base,rel_f)
    remote_f = paste0(remote_base,rel_f)
    print(rel_f)
    if(file.exists(local_f)){
      mtime = file.info(local_f)$mtime
      mtime = stringr::str_remove(mtime,"\\s\\d+:\\d+:\\d+")
      #print(mtime)

      remote_session = check_remote_configuration(auto_connect=TRUE)
      #print(remote_f)
      if(remote_session){
        output = ssh::ssh_exec_internal(ssh_session,paste("stat -L ",remote_f,"| grep Modify"))$stdout

        output = rawToChar(output) %>% stringr::str_extract(.,"\\d+-\\d+-\\d+")
      }else{
        output = file.info(remote_f)$mtime %>% stringr::str_remove("\\s\\d+:\\d+:\\d+")
      }
      remote_time = lubridate::as_date(output)
      local_time = lubridate::as_date(mtime)
      agediff = lubridate::time_length(remote_time - local_time,unit="days")
      if(agediff>0){
        print(paste("Warning! Remote version is",agediff,"days newer than the local file. You probably need to update the following file:"))
        print(rel_f)
      }else{
        message("OK")
      }
    }else{
      if(archive_mode){
        message("local backup of this file doesn't seem to exist:")
        message(local_f)

        copy_no_clobber(remote_f,local_f,force_backup)

      }
    }
  }
}



#' @export
copy_no_clobber = function(from_file,
                           to_file,
                           force = FALSE){

  to_dir = dirname(to_file)
  suppressMessages(suppressWarnings(dir.create(to_dir,recursive = T)))
  print(paste("COPYING",from_file,"TO",to_file))
  if(force){
    file.copy(from_file,to_file)
  }
}



#' @export
get_runs_table = function(seq_type_filter="genome"){
  t_meta = get_gambl_metadata(tissue_status_filter = c("tumour"),seq_type_filter=seq_type_filter) %>%
    dplyr::select(sample_id,patient_id,seq_type,genome_build,pairing_status,unix_group) %>%
    rename("tumour_sample_id"="sample_id")
  n_meta = get_gambl_metadata(tissue_status_filter = c("normal"),seq_type_filter=seq_type_filter) %>%
    dplyr::select(sample_id,patient_id,seq_type,genome_build) %>% rename("normal_sample_id"="sample_id")
  runs_df = left_join(t_meta,n_meta,by=c("patient_id","seq_type","genome_build"))
  #fill in normal_sample_id for unmatched cases
  unmatched_df = get_unmatched_normals(seq_type_filter=seq_type_filter)
  runs_df = left_join(runs_df,unmatched_df,by=c("seq_type","genome_build","unix_group")) %>%
    mutate(normal_sample_id=ifelse(is.na(normal_sample_id.x),normal_sample_id.y,normal_sample_id.x)) %>%
    select(-normal_sample_id.x, -normal_sample_id.y)
  return(runs_df)
}



#' @export
get_template_wildcards = function(parent_key,
                                  template_key){

  if(missing(template_key)){
    wildcard_string = config::get(parent_key)
  }else{
    wildcard_string = config::get(paste0(parent_key,"_wildcards"))[template_key]
  }
  wildcards = stringr::str_split(wildcard_string,",")
  return(unlist(wildcards))
}



#helper function to get the unmatched normals from the main config
#' @export
get_unmatched_normals = function(seq_type_filter){
  a = check_config_value(config::get("unmatched_normal_ids"))
  df = melt(a,value.name="normal_sample_id") %>%
    rename(c("genome_build"="L3","seq_type"="L2","unix_group"="L1")) %>%
    dplyr::filter(seq_type == seq_type_filter)
  return(df)
}



#Helper functions not for export
#' @export
grob_wildcards = function(wildcarded_string){
  wildcards = unlist(stringr::str_extract_all(wildcarded_string,"\\{[^\\{]+\\}"))
  wildcards = stringr::str_remove_all(wildcards,"\\{") %>%  stringr::str_remove_all(.,"\\}")
  return(wildcards)
}
