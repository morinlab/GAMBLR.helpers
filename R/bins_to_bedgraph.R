#' @title Bins to bed graph.
#'
#' @description This function will generate bed file.
#'
#' @param bin_df Data frame with bins.
#' @param min_value Minimum value to retain the bin.
#' @param filename Path to a local file on drive to save the resulting file.
#'
#' @return full path to the file that was written
#'
#' @export
bins_to_bedgraph = function(bin_df,min_value = 3,filename = "test.bed"){
  bed_cols = dplyr::select(bin_df,1,2,3,smoothed_ratio)
  colnames(bed_cols) = c("chr","start","end","value")
  if(any(!grepl("chr",bed_cols[,1]))){
    message("adding chr prefix")
    bed_cols = mutate(bed_cols,chr = paste0("chr",chr))
  }
  bed_cols = dplyr::filter(bed_cols,value > min_value)
  this = dplyr::filter(chr1p,bin_start==27985000)
  bed_cols = mutate(bed_cols,end=format(end, scientific=F),start=format(start, scientific=F))
  write.table(bed_cols,row.names=F,col.names=F,quote=F,sep="\t",file=filename)
  return(filename)
}
