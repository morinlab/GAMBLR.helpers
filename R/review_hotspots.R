
#' @title Annotate curated hotspots
#' 
#' @description Annotate MAF-like data frame with a hot_spot column indicating recurrent mutations.
#'
#' @details This function annotates a MAF-like data frame and will create or overwrite the
#' "hot_spot" column based on curated hotspot rules. Genes for hotspot review are supplied
#' with the `genes_of_interest` parameter.
#' Currently only a few sets of genes are supported, see parameter description for more information and limitations.
#' The desired genome build can be specified with `genome_build` parameter (useful if you loaded the MAF from disk).
#' Obviously, you need to specify the same genome_build as the coordinate system used in the incoming MAF.
#'
#' @param annotated_maf A maf_data object or data frame in MAF format that has hotspots annotated using the function annotate_hotspots().
#' @param genes_of_interest An optional vector of genes for hotspot annotation. By default, the genes present in the curated hotspot tables are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF file. The default is grch37 genome build.
#' @param custom_coordinates A data frame with custom coordinates for the hot spots. 
#' All mutations in any of the regions specified in the data frame will be marked as hot spots.
#' The data frame must have the following columns: "Hugo_Symbol", "chrom", "start", and "end".
#' Optional columns are: "classes" (regex for Variant_Classification) and "alias" (identifier for the hotspot).
#' @param existing_values_action Character. How to handle existing columns.
#' Use "clobber" (default) to always reset "hot_spot" and "mutation_alias"
#' before re-annotating. Use "update" to only fill missing values.
#' @return The same data frame (as given to the `annotated_maf` parameter) with the reviewed columns "hot_spot" and "mutation_alias".
#'
#' @import dplyr
#' @export
#'
#' @examples
#' 
#' coding_capture_maf = GAMBLR.open::get_coding_ssm(
#'   these_samples_metadata = get_gambl_metadata() %>%
#'     dplyr::filter(cohort=="dlbcl_reddy"),
#'   this_seq_type = "capture")
#' ## annotate all supported hotspots
#' coding_capture_maf_anno = annotate_curated_hotspots(coding_capture_maf)
#' 
#' dplyr::filter(coding_capture_maf_anno, hot_spot==TRUE) %>%
#'   dplyr::select(Hugo_Symbol,Start_Position, HGVSp_Short, mutation_alias) %>%
#'   as.data.frame() %>% 
#'   unique()
#' # How does this look for FOXO1?
#' filter(gambl_all_coding,Hugo_Symbol=="FOXO1")  %>% 
#'   dplyr::group_by(mutation_alias) %>% dplyr::count()
#' 



annotate_curated_hotspots = function(annotated_maf,
                           genes_of_interest = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2"),
                           genome_build,
                           custom_coordinates,
                           existing_values_action = "clobber"){
  original_has_maf_class = "maf_data" %in% class(annotated_maf)
  if(missing(genome_build)){
    if("maf_data" %in% class(annotated_maf)){
      genome_build = get_genome_build(annotated_maf)
      #drop our S3 classes because these additional attributes seem to cause some problems when the data is subsequently munged.
      annotated_maf = strip_genomic_classes(annotated_maf)
    }else{
      if("genome_build" %in% colnames(annotated_maf)){
        gb = unique(annotated_maf$genome_build[!is.na(annotated_maf$genome_build)])
        if(length(gb) == 1){
          genome_build = gb
        }else{
          stop("genome_build is required (multiple or missing genome_build values in data)")
        }
      }else if("NCBI_Build" %in% colnames(annotated_maf)){
        gb = unique(annotated_maf$NCBI_Build[!is.na(annotated_maf$NCBI_Build)])
        if(length(gb) == 1){
          genome_build = gb
        }else{
          stop("genome_build is required (multiple or missing NCBI_Build values in data)")
        }
      }else{
        stop("genome_build is required")
      }
    }
  }

  if(!missing(custom_coordinates)){
    #check for the required columns
    if(any(!"Hugo_Symbol" %in% colnames(custom_coordinates))){
      stop("coordinates data frame must have Hugo_Symbol column")
    }
    if(!"chrom" %in% colnames(custom_coordinates)){
      stop("custom_coordinates requires a column named chrom specifying the chromosome of each hot spot region")
    }
    if(!"start" %in% colnames(custom_coordinates) || !"end" %in% colnames(custom_coordinates) ){
      stop("custom_coordinates requires a start and end column specifying the boundaries of each hot spot region")
    }
    coordinates = custom_coordinates
    if(!"classes" %in% colnames(coordinates)){
      coordinates = coordinates %>%
        dplyr::mutate(classes = paste0(vc_nonSynonymous, collapse="|"))
    }
  }else{
    coordinates = get_hotspot_coordinates(genome_build)
  }

  if(!existing_values_action %in% c("clobber", "update")){
    stop("existing_values_action must be one of: clobber, update")
  }

  if(existing_values_action == "clobber"){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(
        hot_spot = "FALSE",
        mutation_alias = paste0(Hugo_Symbol, "_other")
      )
  }else{
    if(!"hot_spot" %in% colnames(annotated_maf)){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = "FALSE")
    }
    if(!"mutation_alias" %in% colnames(annotated_maf) && "hotspot_alias" %in% colnames(annotated_maf)){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(mutation_alias = hotspot_alias)
    }
    if(!"mutation_alias" %in% colnames(annotated_maf)){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(mutation_alias = paste0(Hugo_Symbol, "_other"))
    }else{
      annotated_maf = annotated_maf %>%
        dplyr::mutate(mutation_alias = ifelse(
          is.na(mutation_alias) | mutation_alias == "",
          paste0(Hugo_Symbol, "_other"),
          mutation_alias
        ))
    }
  }
  if("hotspot_alias" %in% colnames(annotated_maf)){
    annotated_maf = annotated_maf %>%
      dplyr::select(-hotspot_alias)
  }
  if(!".row_id" %in% colnames(annotated_maf)){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(.row_id = dplyr::row_number())
  }

  if(!missing(genes_of_interest) && !is.null(genes_of_interest)){
    supported_genes = sort(unique(coordinates$Hugo_Symbol))
    if (length(intersect(supported_genes, genes_of_interest))==0){
        stop(paste0("Currently only ",  paste(supported_genes, collapse=", "),
                    " are supported. Please specify one of these genes."))
    }
    if (length(setdiff(genes_of_interest, supported_genes))>0){
        message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                               " are supported. By default only these genes from the",
                               " supplied list will be reviewed. Reviewing hotspots for genes ",
                               paste(intersect(supported_genes, genes_of_interest),
                                     collapse = ", "), ", it will take a second ...")))
    }
    coordinates = coordinates %>%
      dplyr::filter(Hugo_Symbol %in% genes_of_interest)
  }

  coordinates = coordinates %>%
    dplyr::mutate(
      max = ifelse(start > end, start, end),
      min = ifelse(end < start, end, start)
    ) %>%
    dplyr::mutate(start = min, end = max) %>%
    dplyr::select(-min, -max)

  if("Chromosome" %in% colnames(annotated_maf)){
    annotated_maf = annotated_maf %>%
      dplyr::mutate(Chromosome = as.character(Chromosome))
  }
  if("chrom" %in% colnames(coordinates)){
    coordinates = coordinates %>%
      dplyr::mutate(chrom = as.character(chrom))
  }

  annotated_maf = annotated_maf %>%
    dplyr::filter(Hugo_Symbol %in% coordinates$Hugo_Symbol)

  reviewed_maf = cool_overlaps(
    annotated_maf,
    coordinates,
    columns1 = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position"),
    columns2 = c("Hugo_Symbol", "chrom", "start", "end"),
    type = "any",
    nomatch = TRUE
  )

  original_cols = setdiff(colnames(annotated_maf), ".row_id")

  reviewed_maf = reviewed_maf %>%
    dplyr::mutate(
      .is_hot = !is.na(start) & !is.na(classes) &
        mapply(grepl, pattern = classes, x = Variant_Classification)
    ) %>%
    dplyr::group_by(.row_id) %>%
    dplyr::summarise(
      dplyr::across(all_of(original_cols), ~ dplyr::first(.x)),
      hot_spot = ifelse(any(.is_hot), "TRUE", dplyr::first(hot_spot)),
      mutation_alias = {
        aliases = alias[.is_hot & !is.na(alias)]
        if(length(aliases) > 0) aliases[1] else dplyr::first(mutation_alias)
      },
      .groups = "drop"
    )

  reviewed_maf = reviewed_maf %>%
    dplyr::arrange(.row_id) %>%
    dplyr::select(-.row_id)

  if(original_has_maf_class && exists("create_maf_data", mode = "function")){
    reviewed_maf = create_maf_data(reviewed_maf, genome_build)
  }

  return(reviewed_maf)
}



vc_truncating = c("Nonsense_Mutation","Splice_Site","Frame_Shift_Ins","Frame_Shift_Del")
vc_non_truncating = vc_nonSynonymous[!vc_nonSynonymous %in% vc_truncating]

#' @title Hotspot coordinate definitions
#'
#' @description Return curated hotspot coordinate regions for a given genome build,
#' including per-region variant class patterns and aliases.
#'
#' @param genome_build Reference genome build for the coordinates in the MAF file.
#' Supported values include hg19/grch37/hs37d5/GRCh37 and hg38/grch38/GRCh38.
#' @return A data frame with columns including "Hugo_Symbol", "chrom", "start",
#' "end", "strand", "type", "size", "alias", and "classes".
#'
#' @import dplyr tidyr
#' @export
#'
get_hotspot_coordinates = function(genome_build){
  if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = hotspot_regions_grch37
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = hotspot_regions_hg38
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }
  coordinates = coordinates %>%
    dplyr::rename(Hugo_Symbol = gene)

  #for some reason sometimes this data frame has the higher number as the start.
  coordinates = dplyr::mutate(coordinates,
    max = ifelse(start > end, start, end),
    min = ifelse(end < start, end, start)
  ) %>%
    dplyr::mutate(start = min, end = max) %>%
    dplyr::select(-min, -max)

  coordinates = dplyr::mutate(coordinates, size = end - start + 1)

  # Treat single-position hotspots as "start and above"
  coordinates = coordinates %>%
    mutate(
      end = ifelse(size == 1 & (is.na(strand) | strand != "-"), 300000000, end),
      start = ifelse(size == 1 & strand == "-", 0, start)
    )

  type_to_classes = function(type_value){
    if(is.na(type_value) || type_value == "" || type_value == "all"){
      return(paste0(vc_nonSynonymous, collapse = "|"))
    }
    tokens = unlist(strsplit(type_value, ","))
    tokens = trimws(tokens)
    classes = character(0)
    if("trunc" %in% tokens){
      classes = c(classes, vc_truncating)
    }
    if("missense" %in% tokens){
      classes = c(classes, "Missense_Mutation")
    }
    if("inframe" %in% tokens){
      classes = c(classes, "In_Frame_Del", "In_Frame_Ins")
    }
    if("stoploss" %in% tokens){
      classes = c(classes, "Nonstop_Mutation")
    }
    if("startloss" %in% tokens){
      classes = c(classes, "Translation_Start_Site")
    }
    if(length(classes) == 0){
      classes = vc_nonSynonymous
    }
    paste0(unique(classes), collapse = "|")
  }

  coordinates = coordinates %>%
    dplyr::mutate(classes = vapply(type, type_to_classes, character(1)))
  
  return(coordinates)
  
}

#' @title Review Hotspots.
#'
#' @description Annotate MAF-like data frome with a hot_spot column indicating
#'      recurrent mutations.
#'
#' @details This function takes an annotated MAF (the output of
#'      annotate_hotspots) and reviews an existing column, "hot_spot", in the
#'      provided data frame. Genes for hotspot review are supplied with the
#'      `genes_of_interest` argument. Currently only a few sets of genes are
#'      supported, see parameter description for more information. The desired
#'      genome build can be specified with `genome_build` argument and is
#'      expected to be the same as the incoming MAF. The review of hotspots is
#'      gene-specific and manually builds on top of the provided oncodrive
#'      annotations. Specifically, regardless of the oncodrive annotations,for
#'      FOXO1 all mutations affecting "p.M1?" will be annotated as hotspots. For
#'      MYD88, all mutations affecting "p.L273P" and "p.L265P" will be annotated
#'      as hotspots. For CREBBP, all missense mutations in the KAT domain will
#'      be annotated hotspots. For NOTCH1 and NOTCH2, all PEST domain mutations
#'      will become annotated hotspots. For CD79B, the truncating mutations will
#'      become annotated hotspots. Finally, for EZH2 all mutations within the
#'      7:148508764-148506238 (grch37, also supported in hg38 projection) will
#'      become annotated as hotspots.
#'
#' @param annotated_maf A data frame in MAF format that has hotspots annotated
#'      using the function annotate_hotspots().
#' @param genes_of_interest A vector of genes for hotspot review. Currently only
#'      FOXO1, MYD88, CREBBP, NOTCH1, NOTCH2, CD79B and EZH2 are supported.
#' @param genome_build Reference genome build for the coordinates in the MAF
#'      file. The default is grch37 genome build.
#' @param custom_coordinates A data frame with custom coordinates for the hot
#'      spots. All mutations in any of the regions specified in the data frame
#'      will be marked as hot spots. The data frame must have the following
#'      columns: "Hugo_Symbol", "chrom", "start", and "end".
#'
#' @return data frame
#'
#' @import dplyr
#' @export
#'
#' @examples
#' hot_ssms = review_hotspots(annotate_hotspots(get_coding_ssm(this_seq_type = "genome")),
#'                            genes_of_interest = c("CREBBP"))
#'
review_hotspots = function(annotated_maf,
                           genes_of_interest = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2"),
                           genome_build,
                           custom_coordinates){
  if(missing(genome_build)){
    if("maf_data" %in% class(annotated_maf)){
      genome_build = get_genome_build(annotated_maf)
      #drop our S3 classes because these additional attributes seem to cause some problems when the data is subsequently munged.
      annotated_maf = strip_genomic_classes(annotated_maf)
    }else{
      stop("genome_build is required")
    }
  }

  # define the list of genes currently supported for review
  supported_genes = c("FOXO1", "MYD88", "CREBBP", "NOTCH1", "NOTCH2", "CD79B", "EZH2")

  if(!missing(custom_coordinates)){
    #check for the required columns
    if(any(!"Hugo_Symbol" %in% colnames(custom_coordinates))){
      stop("coordinates data frame must have Hugo_Symbol column")
    }
    if(!"chrom" %in% colnames(custom_coordinates)){
      stop("custom_coordinates requires a column named chrom specifying the chromosome of each hot spot region")
    }
    if(!"start" %in% colnames(custom_coordinates) || !"end" %in% colnames(custom_coordinates) ){
      stop("custom_coordinates requires a start and end column specifying the boundaries of each hot spot region")
    }
    coordinates = custom_coordinates
  }else if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = hotspot_regions_grch37 %>%
      dplyr::rename(Hugo_Symbol = gene)
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = hotspot_regions_hg38 %>%
      dplyr::rename(Hugo_Symbol = gene)
  }else{
    stop("The genome build specified is not currently supported. Please provide MAF file in one of the following cordinates: hg19, grch37, hs37d5, GRCh37, hg38, grch38, or GRCh38")
  }
  
  if(!missing(custom_coordinates)){
    if(!missing(genes_of_interest)){
      custom_coordinates = filter(custom_coordinates,Hugo_Symbol %in% genes_of_interest)
    }
    print("HERE!")
    custom_coordinates = custom_coordinates %>% dplyr::select(-Hugo_Symbol)
    #everything is just naive coordinate range-based
    annotated_maf = cool_overlaps(annotated_maf,
                                  custom_coordinates,
                                  columns1 = c("Chromosome","Start_Position","End_Position"),
                                  columns2 = c("chrom","start","end"),
                                  type = "any",
                                  nomatch = TRUE) 
      print(table(is.na(annotated_maf$start)))
      annotated_maf = annotated_maf %>% 
        mutate(hot_spot = ifelse(is.na(start),FALSE,TRUE))
      message("Hotspot review complete")
      return(annotated_maf)
  }
  # check that at least one of the currently supported genes are present
  if (length(intersect(supported_genes, genes_of_interest))==0){
      stop(paste0("Currently only ",  paste(supported_genes, collapse=", "), 
                  " are supported. Please specify one of these genes."))
  }
  # notify user that there is limited number of genes currently supported
  

  if (length(setdiff(genes_of_interest, supported_genes))>0){
      message(strwrap(paste0("Currently only ", paste(supported_genes, collapse=", "),
                             " are supported. By default only these genes from the",
                             "supplied list will be reviewed. Reviewing hotspots for genes ",
                             paste(intersect(supported_genes, genes_of_interest),
                                   collapse = ", "), ", it will take a second ...")))
  }

  if("FOXO1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "FOXO1" &
                                        HGVSp_Short == "p.M1?",
                                        "TRUE", hot_spot))
  }
  get_coordinate = function(gene, type_pattern = NULL){
    coord = coordinates %>% dplyr::filter(Hugo_Symbol == gene)
    if(!is.null(type_pattern)){
      coord = coord %>% dplyr::filter(grepl(type_pattern, type))
    }
    if(nrow(coord) == 0){
      stop(paste0("No hotspot coordinates found for ", gene))
    }
    coord[1, , drop = FALSE]
  }


  if("CREBBP" %in% genes_of_interest){
      crebbp = get_coordinate("CREBBP", "missense")
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CREBBP" &
                                        Start_Position > crebbp$start &
                                        End_Position < crebbp$end &
                                        Variant_Classification == "Missense_Mutation",
                                        "TRUE", hot_spot))
  }
  if("EZH2" %in% genes_of_interest){
      ezh2 = get_coordinate("EZH2")
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "EZH2" &
                                        Start_Position > ezh2$start &
                                        End_Position < ezh2$end,
                                        "TRUE", hot_spot))
  }
  if("MYD88" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "MYD88" &
                                        HGVSp_Short %in% c("p.L273P", "p.L265P"),
                                        "TRUE", hot_spot))
  }
  if("NOTCH1" %in% genes_of_interest){
      notch1 = get_coordinate("NOTCH1")
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH1" &
                                        Start_Position < notch1$start,
                                        "TRUE", hot_spot))
  }
  if("NOTCH2" %in% genes_of_interest){
      notch2 = get_coordinate("NOTCH2")
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH2" &
                                        Start_Position < notch2$start,
                                        "TRUE", hot_spot))
  }

  if("CD79B" %in% genes_of_interest){
      cd79b_trunc = get_coordinate("CD79B", "trunc")
      cd79b_nontrunc = get_coordinate("CD79B", "missense|inframe")
      truncating_variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site")
      annotated_maf = annotated_maf %>%
         dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                         Start_Position < cd79b_trunc$start &
                                         Variant_Classification %in% truncating_variants,
                                         "TRUE", hot_spot)) %>%
          dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                          Start_Position < cd79b_nontrunc$start &
                                          ! Variant_Classification %in% truncating_variants,
                                          "TRUE", hot_spot))
  }
  return(annotated_maf)
}
