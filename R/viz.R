# global set of aliases for finding specific sets of colours
#' @export
colour_aliases = list("COO_consensus" = "coo", "COO" = "coo", "DHITsig_consensus" = "coo",
                      "pathology" = "pathology", "analysis_cohort" = "pathology", "group" = "pathology",
                      "FL_group" = "FL", "lymphgen" = "lymphgen", "lymphgen_with_cnv" = "lymphgen",
                      "bcl2_ba" = "pos_neg", "BCL2_status" = "pos_neg", "myc_ba" = "pos_neg",
                      "bcl6_ba" = "pos_neg", "EBV_status_inf"="EBV_status",
                      "manta_BCL2_sv" = "pos_neg", "manual_BCL2_sv" = "pos_neg", "manta_MYC_sv" = "pos_neg")



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



#' INTERNAL FUNCTION called by [GAMBLR::plot_sample_circos] and [GAMBLR::get_mutation_frequency_bin_matrix], not meant for out-of-package usage.
#' Assign a colour palette to metadata columns automatically and consistently.
#'
#' @param metadataColumns Names of the metadata columns to assign colours for.
#' @param these_samples_metadata Metadata for just the samples you need colours for.
#' @param column_alias A list of column_names with values indicating what gambl colour scheme they should use (e.g. pos_neg, pathology, lymphgen).
#' @param as_vector Boolean statement that is set to TRUE per default.
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param annoAlpha Optional alpha to apply to annotation colours.
#'
#' @return Either a vector or list of colours.
#'
#' @import dplyr ggsci
#'
#' @noRd
#'
#' @examples
#' #get metadata
#' all_meta = get_gambl_metadata()
#'
#' #get colours
#' all_cols = map_metadata_to_colours(metadataColumns = c("lymphgen",
#'                                                        "pathology",
#'                                                        "genetic_subgroup"),
#'                                    these_samples_metadata = all_meta,
#'                                    column_alias = list("nothing" = "FL"),
#'                                    as_vector = F)
#'
#' @export
map_metadata_to_colours = function(metadataColumns,
                                   these_samples_metadata,
                                   column_alias = list(),
                                   as_vector = TRUE,
                                   verbose = FALSE,
                                   annoAlpha = 1){

  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours = get_gambl_colours()
  colours = list()
  colvec = c()


  aliases = c(colour_aliases, column_alias)
  for(column in metadataColumns){
    this_value = these_samples_metadata[[column]]
  options = this_value
    if(verbose){
      print(">>>>>>>")
      message("finding colour for", this_value)
      print("<<<<<<<")
    }
    if(column %in% names(aliases)){
      key = aliases[[column]]
      if(verbose){
        print(paste("using alias to look up colours for", column))
        message(paste("using", key, "for", column))
      }
      these = get_gambl_colours(classification = key)
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      if(verbose){
        message("adding:", these[this_value])
      }
    }else if(column == "sex"){
      these = get_gambl_colours("sex", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      message("adding:", these[this_value])
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){

      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "in clinical"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if(("positive" %in% options | "POS" %in% options) & length(options)<4){
      if(verbose){
        print("using pos_neg")
      }

      these = get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]

      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if("GCB" %in% options){
      these = get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec,these)
    }else if(column %in% c("pathology")){
      these = get_gambl_colours(column, alpha = annoAlpha)

      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(column == "HMRN"){
      these = get_gambl_colours("hmrn", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(length(levels(options)) > 15){

      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)

      colours[[column]] = these
      colvec = c(colvec, these)
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }
  }
  if(as_vector){
    return(colvec)
  }
  return(colours)
}



#' This function doesn't do anything yet, thus, it's not currently exported to NAMESPACE.
#'
#' @param mafs TODO
#' @param this_sample_id TODO
#' @param genes TODO
#' @param show_noncoding TODO
#' @param detail TODO
#'
#' @noRd
#'
#' @export
plot_multi_timepoint = function(mafs,
                                this_sample_id,
                                genes,
                                show_noncoding = FALSE,
                                detail){

  tp = c("A","B","C")
  title = paste(sample_id, detail, sep = "\n")
  i = 1
  for (i in c(1:length(mafs))){
    maf_file = mafs[i]
    time_point = tp[i]
    print(paste(maf_file, time_point))
  }
  if(length(mafs)==2){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])

    A.maf = A.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count)) %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 1) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    B.maf = B.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count)) %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 2) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    all.maf = rbind(A.maf, B.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf, shm_regions)
    }
    else{
      coding.maf = dplyr::filter(all.maf, !Variant_Classification %in% c("Silent", "RNA", "IGR", "Intron", "5'Flank", "3'Flank", "5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf, (Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }
    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf, coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf, coord %in% B.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords), "category"] = "not-A"
    coding.maf[which(coord %in% B.zero.coords), "category"] = "not-B"

    #actually this is changed in eitehr direction, not just gained
    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" )

    #just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point !=2) %>%
    #  mutate(time_point = time_point +0.4)

    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point == 2 & VAF > 0) %>%
      mutate(time_point = time_point +0.4)

    print(just_gained_lg)

    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)

    ggplot(coding.maf, aes(x = time_point, y = VAF, group = coord, colour = category)) +
      geom_point() +
      geom_line(alpha = 0.5) +
      geom_text_repel(data = just_gained_lg, aes(label = Hugo_Symbol), size = 4,segment.linetype = 0) +
      geom_text_repel(data = just_trunk, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  if(length(mafs)==3){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])
    C.maf = fread_maf(mafs[3])

    A.maf = A.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 1) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    B.maf = B.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 2) %>%
      mutate(coord = paste(Chromosome, Start_Position,sep = ":"))

    C.maf = C.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 3) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    all.maf = rbind(A.maf, B.maf, C.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf, shm_regions)
    }
    else{
      coding.maf = dplyr::filter(all.maf, !Variant_Classification %in% c("Silent", "RNA", "IGR", "Intron", "5'Flank", "3'Flank", "5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf, (Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }

    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf, coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf, coord %in% B.zero.coords)

    C.rows = which(coding.maf$time_point==3)
    C.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==3 & VAF == 0), "coord"]))
    C.zero = dplyr::filter(coding.maf, coord %in% C.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords), "category"] = "not-A"
    coding.maf[which(coord %in% B.zero.coords), "category"] = "not-B"
    coding.maf[which(coord %in% C.zero.coords), "category"] = "not-C"

    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" )

    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" & time_point ==3) %>%
      mutate(time_point = time_point +0.4)

    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)

    ggplot(coding.maf, aes(x = time_point, y = VAF, group = coord, colour = category)) +
      geom_point() +
      geom_line(alpha = 0.3) +
      geom_text_repel(data = just_gained_lg, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      geom_text_repel(data = just_trunk, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  return(FALSE)
}



#' @title VAF Heatmap
#'
#' @description Plot a heatmap comparing the VAF of mutations in T1/T2 pairs.
#'
#' @details Currently unfinished plotting function. Thus, I have removed it from export until it's in a state where it can be included in GAMBLR.
#' Parameter descriptions need to be updated so that the origin of the incoming data is clear.
#' Examples would also need to be added before this function gets exported into NAMESPACE.
#'
#' @param maf1 Data frame of simple somatic mutations at time point A.
#' @param maf2 Data frame of simple somatic mutations at time point B.
#' @param vafcolname Name of variable that holds VAF in maf. If not supplied, vaf will be calcualted.
#' @param patient_colname Column name that holds patient name (default is "patient_id").
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata).
#' @param sortByColumns Which of the metadata to sort on for the heatmap.
#' @param metadata_columns A vector containing the categorical column names you want to plot below.
#' @param gene_orientation Where genes would be plotted. Default is "bottom".
#' @param annotate_zero Indicate a variant that had VAF = 0 in one of the two time points. Default is FALSE.
#' @param genes An optional vector of genes to restrict your plot to.
#' @param top_n_genes How many genes to be added to the plot.
#' @param drop_unless_lowvaf Will drop some genes where VAF is low, default is FALSE.
#' @param vaf_cutoff_to_drop Which VAF cut-off value to use when dropping variants before plotting.
#' @param cluster_columns Boolean statement for clustering by columns, defaults to FALSE.
#' @param cluster_rows Boolean statement for clustering by rows, defaults to FALSE.
#'
#' @return Nothing
#' @noRd
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr circlize ComplexHeatmap tibble
#'
#' @export
plot_mutation_dynamics_heatmap = function(maf1,
                                          maf2,
                                          vafcolname,
                                          patient_colname = "patient_id",
                                          these_samples_metadata,
                                          sortByColumns,
                                          metadata_columns = c("sample_id"),
                                          gene_orientation = "bottom",
                                          annotate_zero = FALSE,
                                          genes,
                                          top_n_genes,
                                          drop_unless_lowvaf = FALSE,
                                          vaf_cutoff_to_drop = 0.05,
                                          cluster_columns = FALSE,
                                          cluster_rows = FALSE){

  if(missing(vafcolname)){
    t1_pair_mafdat = mutate(maf1, vaf = t_alt_count / (t_alt_count + t_ref_count))
    t2_pair_mafdat = mutate(maf2, vaf = t_alt_count / (t_alt_count + t_ref_count))
  }else{
    t1_pair_mafdat[,"vaf"] = t1_pair_mafdat[,vafcolname]
    t2_pair_mafdat[,"vaf"] = t2_pair_mafdat[,vafcolname]
  }
  if(!missing(sortByColumns)){
    these_samples_metadata = arrange(these_samples_metadata, across(sortByColumns))
  }

  t1_pair_mafdat = t1_pair_mafdat %>%
    dplyr::rename("patient_id" = patient_colname)

  t2_pair_mafdat = t2_pair_mafdat %>%
    dplyr::rename("patient_id" = patient_colname)

  these_samples_metadata = these_samples_metadata %>%
    dplyr::rename("patient_id" = patient_colname)

  both_vaf_all = full_join(t1_pair_mafdat, t2_pair_mafdat, by = c("patient_id", "Start_Position")) %>%
    dplyr::select(patient_id, HGVSp_Short.x, Hugo_Symbol.x, Hugo_Symbol.y, Tumor_Sample_Barcode.x,
                  Tumor_Sample_Barcode.y, HGVSp_Short.y, vaf.x, vaf.y, hot_spot.x, hot_spot.y) %>%
    mutate(ANNO = ifelse(is.na(vaf.x), HGVSp_Short.y, HGVSp_Short.x)) %>%
    mutate(GENE = ifelse(is.na(vaf.x), Hugo_Symbol.y, Hugo_Symbol.x)) %>%
    mutate(VAF1 = ifelse(is.na(vaf.x), 0, vaf.x)) %>%
    mutate(VAF2 = ifelse(is.na(vaf.y), 0, vaf.y)) %>%
    mutate(hot_spot.x = ifelse(is.na(hot_spot.x), 0, 1)) %>%
    mutate(hot_spot.y = ifelse(is.na(hot_spot.y), 0, 1)) %>%
    mutate(HOTSPOT = ifelse(hot_spot.x == 1 | hot_spot.y == 1, 1, 0)) %>%
    dplyr::filter(ANNO != "") %>%
    mutate(Mutation = paste(GENE, ANNO, sep = "_")) %>%
    dplyr::select(patient_id, GENE, Mutation, VAF1, VAF2, ANNO, HOTSPOT)

  if(drop_unless_lowvaf){
    both_vaf_all = dplyr::filter(both_vaf_all, VAF1 < vaf_cutoff_to_drop | VAF2 < vaf_cutoff_to_drop)
  }
  if(!missing(genes)){
    both_vaf_all = dplyr::filter(both_vaf_all, GENE %in% genes)
  }
  both_vaf_all = mutate(both_vaf_all, unique_id = paste(patient_id, Mutation, sep = "_")) %>%
    mutate(fold_change = log(VAF2 + 0.1)-log(VAF1 + 0.1)) %>%
    group_by(patient_id, GENE) %>%
    mutate(Number = paste(GENE, row_number(), sep = "_"))

  print(head(both_vaf_all))
  h = both_vaf_all %>%
    select(patient_id, Number, fold_change) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = fold_change) %>%
    column_to_rownames("patient_id")

  both_vaf_all = mutate(both_vaf_all, minVAF = ifelse(VAF1 < VAF2, VAF1, VAF2))
  zeroes = both_vaf_all %>%
    select(patient_id, Number, minVAF) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = minVAF) %>%
    column_to_rownames("patient_id")

  hotspots = both_vaf_all %>%
    select(patient_id, Number, HOTSPOT) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = HOTSPOT) %>%
    column_to_rownames("patient_id")

  zeroes[is.na(zeroes)] = 0.001
  hotspots[is.na(hotspots)] = -1
  hotspots[hotspots == 0] = -1
  hotspots[hotspots == 1] = 0
  zeroes[zeroes > 0] = 1
  print(head(hotspots))

  these_samples_metadata_rn = dplyr::filter(these_samples_metadata, patient_id %in% rownames(h)) %>%
    select(all_of(c("patient_id", metadata_columns))) %>%
    column_to_rownames("patient_id")

  la = HeatmapAnnotation(df = as.data.frame(these_samples_metadata_rn), which = "row")
  ta = HeatmapAnnotation(df = as.data.frame(these_samples_metadata_rn), which = "column")

  h[is.na(h)] = 0
  cs = colSums(zeroes)
  ordered = names(cs[order(cs)])
  if(!missing(top_n_genes)){
    print(head(ordered, top_n_genes))
    genes = ordered[c(1:top_n_genes)]
  }
  col_fun = colorRamp2(c(0, 1), c("white", "red"))

  if(gene_orientation == "bottom"){
    if(annotate_zero){

    }else{
      H = Heatmap(h[rownames(these_samples_metadata_rn),], cluster_rows = F, cluster_columns = F, left_annotation = la)
    }
  }else{
    if(!missing(top_n_genes)){
      these_zeroes = t(zeroes[rownames(these_samples_metadata_rn), genes])
      these_zeroes = t(hotspots[rownames(these_samples_metadata_rn), genes])
      to_show = t(h[rownames(these_samples_metadata_rn), genes])

    }else{
      these_zeroes = t(zeroes[rownames(these_samples_metadata_rn),])
      these_zeroes = t(hotspots[rownames(these_samples_metadata_rn),])
      to_show = t(h[rownames(these_samples_metadata_rn),])
    }
    if(annotate_zero){
      H = Heatmap(to_show, layer_fun = function(j, i, x, y, width, height, fill) {
                                                v = pindex(these_zeroes, i, j)
                                                l = v == 0
                                                grid.points(x[l], y[l], pch = 16, size = unit(1, "mm"), gp = gpar(col = "white"))
    },
      cluster_columns = cluster_columns, cluster_rows = cluster_rows, bottom_annotation = ta)
      print("HERE")
    }else{
      H = Heatmap(t(h[rownames(these_samples_metadata_rn),]), cluster_rows = cluster_rows, cluster_columns = cluster_columns, bottom_annotation = ta)
    }
  }
  return(H)
}



#' @title Trim Scale Expressions.
#'
#' @description INTERNAL HELPER FUNCTION called by prettyOncoplot, not meant for out-of-package usage.
#'
#' @details INTERNAL FUNCTION called by prettyOncoplot, not meant for out-of-package usage.
#'
#' @param x Numeric value (of expression) to be trimmed.
#'
#' @return Numeric value.
#'
#' @noRd
#'
#' @examples
#' trimmed = trim_scale_expression(2)
#'
#' @export
trim_scale_expression = function(x){
  quants = unname(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  x = ifelse(x < quants[1], quants[1], x)
  x = ifelse(x > quants[2], quants[2], x)
  x = (x - quants[1]) / (quants[2] - quants[1])
  return(x)
}
