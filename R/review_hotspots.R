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

  # check genome build because CREBBP coordinates are hg19-based or hg38-based
  if(!missing(custom_coordinates)){
    #check for the required columns
    if(any(!"Hugo_Symbol" %in% colnames(custom_coordinates))){
      stop("coordinates data frame must have Hugo_Symbol column")
    }
    if(!"chrom" %in% colnames(custom_coordinates)){
      stop("custom_coordinates requires a chrom column specifying the chromosome of each hot spot region")
    }
    if(!"start" %in% colnames(custom_coordinates) || !"end" %in% colnames(custom_coordinates) ){
      stop("custom_coordinates requires a start and end column specifying the boundaries of each hot spot region")
    }
    coordinates = custom_coordinates
  }else if (genome_build %in% c("hg19", "grch37", "hs37d5", "GRCh37")){
    coordinates = GAMBLR.data::hotspot_regions_grch37 %>% rownames_to_column("Hugo_Symbol")
  }else if(genome_build %in% c("hg38", "grch38", "GRCh38")){
    coordinates = GAMBLR.data::hotspot_regions_hg38 %>% rownames_to_column("Hugo_Symbol")
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

  if("CREBBP" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CREBBP" &
                                        Start_Position > coordinates["CREBBP", "start"] &
                                        End_Position < coordinates["CREBBP", "end"] &
                                        Variant_Classification == "Missense_Mutation",
                                        "TRUE", hot_spot))
  }
  if("EZH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "EZH2" &
                                        Start_Position > coordinates["EZH2", "start"] &
                                        End_Position < coordinates["EZH2", "end"],
                                        "TRUE", hot_spot))
  }
  if("MYD88" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "MYD88" &
                                        HGVSp_Short %in% c("p.L273P", "p.L265P"),
                                        "TRUE", hot_spot))
  }
  if("NOTCH1" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH1" &
                                        Start_Position < coordinates["NOTCH1", "start"],
                                        "TRUE", hot_spot))
  }
  if("NOTCH2" %in% genes_of_interest){
      annotated_maf = annotated_maf %>%
        dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "NOTCH2" &
                                        Start_Position < coordinates["NOTCH2", "start"],
                                        "TRUE", hot_spot))
  }

  if("CD79B" %in% genes_of_interest){
      truncating_variants = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Region", "Splice_Site")
      annotated_maf = annotated_maf %>%
         dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                         Start_Position < coordinates["CD79B_trunc", "start"] &
                                         Variant_Classification %in% truncating_variants,
                                         "TRUE", hot_spot)) %>%
          dplyr::mutate(hot_spot = ifelse(Hugo_Symbol == "CD79B" &
                                          Start_Position < coordinates["CD79B_NONtrunc", "start"] &
                                          ! Variant_Classification %in% truncating_variants,
                                          "TRUE", hot_spot))
  }
  return(annotated_maf)
}
