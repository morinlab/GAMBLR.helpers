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
