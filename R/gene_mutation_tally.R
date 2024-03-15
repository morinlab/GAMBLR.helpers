#' @title Tally coding mutations load.
#'
#' @description Calculate the number of coding mutations per grouping. Silent
#' variants are excluded from this calculation. The incoming maf may contain
#' non-coding variants, which will be excluded from tallying.
#'
#' @param maf_df Data frame of simple somatic mutations in maf format. Required
#'      parameter.
#' @param these_samples_metadata Data frame with metadata. Must contain sample
#'      identifiers in the `sample_id` column and column that will be used
#'      to calulate the mutation load. All other columns are ignored. Required
#'      parameter.
#' @param these_genes Vector of hugo symbols for genes to be considered for the
#'      tallying of mutations. Required parameter.
#' @param grouping_variable Column in the metadata that will be used as grouping
#'      variable. By default, the `cohort` is used.
#'
#' @return data frame
#'
#' @export
gene_mutation_tally = function(maf_df,these_samples_metadata,these_genes,grouping_variable="cohort"){
  meta = dplyr::select(these_samples_metadata,sample_id,{{grouping_variable}})
  maf_filt = dplyr::filter(maf_df,Hugo_Symbol %in% these_genes, Variant_Classification %in% coding_class) %>%
    dplyr::filter(Variant_Classification !="Silent")
  meta_anno = left_join(maf_filt,meta,by=c("Tumor_Sample_Barcode"="sample_id")) %>%
    group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head()
  meta_anno = left_join(maf_filt,meta,by=c("Tumor_Sample_Barcode"="sample_id")) %>%
    group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head() %>%
    ungroup()
  denom = meta %>% group_by(!!sym(grouping_variable)) %>% tally() %>% dplyr::rename(c("total"="n"))
  meta_anno_tally = group_by(meta_anno,Hugo_Symbol,!!sym(grouping_variable)) %>% tally()
  meta_anno_tally = full_join(meta_anno_tally,denom) %>% mutate(frequency=100*n/total)
  return(meta_anno_tally)
}
