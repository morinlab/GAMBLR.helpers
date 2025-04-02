#' Add reciprocal events to a bedpe file
#'
#' @param bed_df A bed or bedpe data frame (e.g. containing SVs from SVAR or Manta)
#' @param shift_by The positions of START_A and END_A will be increased by this value (default 30). When combining this with decorate_bed, this offset ensures IGV shows the intended colour instead of the colour of the reciprocal event.
#' 
#' @return data frame with each row duplicated and the contents of CHROM_A, START_A etc swapped with CHROM_B, START_B
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' recip_bed = reciprocate_bedpe(in_bed)
#' coloured_bed = decorate_bed(recip_bed,
#'                             colour_mapping = get_gambl_colours("chromosome",
#'                                                                 as_rgb_string = T))
#'}
reciprocate_bedpe = function(bed_df,shift_by=30){
  swapped_bed = bed_df
  swapped_bed[,c(1:3)]=bed_df[,c(4:6)]
  swapped_bed[,c(4:6)]=bed_df[,c(1:3)]
  colnames(swapped_bed) = colnames(bed_df)
  swapped_bed = mutate(swapped_bed,
                       START_A = START_A + shift_by,
                       END_A = END_A + shift_by)
  merged_bed = bind_rows(swapped_bed,bed_df)
  return(merged_bed)
}


#' Decorate a bed or bedpe file
#'
#' @param bed_df A bed or bedpe data frame
#' @param colour_by Specify the name of the column used to define the colours
#' @param colour_mapping Named vector specifying the colour for each possible value in the colour_by column 
#' @param new_column_name The name to use for the colour column (default is 'color' for IGV compatability)
#'
#' @return data frame with one additional column containing colours for each row
#' @export
#'
#' @examples
#' \dontrun{
#' in_bed = get_combined_sv(projection="hg38")
#' # TODO: update decorate_bed to handle non-prefixed chromosome names
#' }
#' 
#' in_bed = get_combined_sv(projection="hg38")
#' coloured_bed = decorate_bed(in_bed,
#'                             colour_mapping = get_gambl_colours("chromosome",
#'                                                               as_rgb_string = T))
decorate_bed = function(bed_df,
                        colour_by="CHROM_B",
                        colour_mapping,
                        new_column_name="color"){
  unique_vals = pull(bed_df,colour_by) %>% unique()
  if(!any(unique_vals %in% names(colour_mapping))){
    print(unique_vals)
    stop("none of the names in your colour_mapping vector match the values in the specified colour_by column")
  }
  #create a data frame for each colour/name pair
  col_df = data.frame(colour_mapping) %>% 
    rownames_to_column(colour_by) %>%
    dplyr::rename((!!new_column_name) := "colour_mapping")
  decorated_bed = left_join(bed_df,col_df,by=colour_by)
  return(decorated_bed)
}


