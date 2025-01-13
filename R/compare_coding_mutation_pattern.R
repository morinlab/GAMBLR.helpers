#' @title Compare pattern of coding mutations.
#'
#' @description TODO.
#'
#' @param maf_df1 TODO.
#' @param maf_df2 TODO.
#' @param gene TODO.
#'
#' @return list
#'
#' @import dplyr
#' @export
compare_coding_mutation_pattern = function(maf_df1,maf_df2,gene){
  if(missing(maf_df1) | missing(maf_df2)){
    stop("must provide two data frames containing mutations you would like to compare")
  }
  if(missing(gene)){
    stop("Must provide the Hugo_Symbol of a single gene that is present in both maf files")
  }
  missense_positions1 = dplyr::filter(maf_df1,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% gsub("p\\.\\w", "", .) %>% regmatches(., regexpr("\\d+", .)) %>% as.numeric()
  missense_positions2 = dplyr::filter(maf_df2,Hugo_Symbol==gene,!Variant_Classification %in% c("Silent","Splice_Site","Splice_Region"),Variant_Type=="SNP") %>%
    pull(HGVSp_Short) %>% gsub("p\\.\\w", "", .) %>% regmatches(., regexpr("\\d+", .)) %>% as.numeric()
 if(length(missense_positions1)==0 | length(missense_positions2)==0 ){
   message(paste("no mutations for",gene,"in one or both data sets"))
   return(list(kl=15))
 }
  #generate range of amino acids based on what we can infer from the MAF (not ideal)
  max_pos = max(c(missense_positions1,missense_positions2))
  full_df = data.frame(position=c(1:max_pos))
  df1 = data.frame(position=missense_positions1) %>% group_by(position) %>% tally() %>% rename("group1"="n")
  df2 = data.frame(position=missense_positions2) %>% group_by(position) %>% tally() %>% rename("group2"="n")
  full_df = left_join(full_df,df1) %>% mutate(group1=ifelse(is.na(group1),0,group1))
  full_df = left_join(full_df,df2) %>% mutate(group2=ifelse(is.na(group2),0,group2))
  # convert to the format needed by KL
  all_counts = dplyr::select(full_df,-position) %>% t()
  all_counts[1,]=all_counts[1,]/sum(all_counts[1,])
  all_counts[2,]=all_counts[2,]/sum(all_counts[2,])

  # Normalize the rows to turn counts into probabilities
  P <- all_counts[1, ] / sum(all_counts[1, ])
  Q <- all_counts[2, ] / sum(all_counts[2, ])

  kl_out <- kl_divergence(P, Q)

  return(list(df=full_df,kl=unname(kl_out)))
}
