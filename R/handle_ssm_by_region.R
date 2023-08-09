#' @title Handle SSM by Region
#'
#' @description Subset an incoming MAF to desired regions.
#'
#' @details This helper function is used to subset an incoming MAF to specific region.
#' This function was developed to accommodate non-GSC-access users to use plotting functions that, 
#' previously internally called GABLR core functions (GSC access required).
#' Thsi function takes an already loaded MAF (`this_maf`) and subset this based on the regions provided,
#' either in region format (chrX:1234-5678) withh the `region` parameter.
#' Or, the suer can seperatly specify the regions of interest using the following parameters;
#' `chromosome`, `qstart`, and, `qend`.
#'
#' @param this_maf An already loaded MAF or MAF-like object to subset to regions of interest. This is a required paremter.
#' @param chromosome The chromosome you are restricting to (with or without a chr prefix).
#' @param qstart Query start coordinate of the range you are restricting to.
#' @param qend Query end coordinate of the range you are restricting to.
#' @param region Region formatted like chrX:1234-5678 instead of specifying chromosome, start and end separately.
#'
#' @return A MAF that has been subset to the regions specified.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' my_maf = GAMBLR.data::sample_data$grch37$maf
#' my_mutations = handle_ssm_by_region(this_maf = my_maf, region = "chr8:128,723,128-128,774,067")
#'
handle_ssm_by_region = function(this_maf,
                                chromosome,
                                qstart,
                                qend,
                                region = ""){
  
  if(!region == ""){
     region = gsub(",", "", region)
     split_chunks = unlist(strsplit(region, ":"))
     chromosome = split_chunks[1]
     startend = unlist(strsplit(split_chunks[2], "-"))
     qstart = as.numeric(startend[1])
     qend = as.numeric(startend[2])
  }
  
  #ensure that the region is specified according to the projection of incoming maf (chr prefix or not)
  if(all(str_detect(this_maf$Chromosome, "chr"))){
    if(!str_detect(chromosome, "chr")){
      chromosome = paste0("chr", chromosome) %>% as.character()
    }
  }else{
    if(str_detect(chromosome, "chr")){
      chromosome = gsub("chr", "", chromosome)
    }
  }

  #subset MAF to region
  muts_region = dplyr::filter(this_maf, Chromosome == chromosome & Start_Position > qstart & Start_Position < qend)
  
  return(muts_region)
}
