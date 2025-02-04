#' @title Get GAMBL Colours.
#'
#' @description Get GAMBL colour schemes for annotating figures.
#'
#' @details This function was designed to retrieve specified GAMBL colour palettes.
#' By default, this function returns all the colours currently available.
#' The user can easily specify what classification to return colors for with the `classification` parameter.
#' It is also possible to return any given colour in different formats.
#' To do so, refer to the Boolean arguments; `as_list` and `as_dataframe`.
#' For more information regarding the available colours, refer to the utilities vignette.
#'
#' @param classification Optionally request only colours for pathology, lymphgen, mutation or copy_number.
#' @param alpha Alpha of plotted colours.
#' @param as_list Boolean parameter controlling the format of the return. Default is FALSE.
#' @param as_dataframe Boolean parameter controlling the format of the return. Default is FALSE.
#' @param return_available Set to TRUE for returning all available colours. Default is FALSE.
#' @param verbose Default is FALSE
#' @param as_rgb_string Set to TRUE if you want RGB triples instead of hex codes (e.g. "#555FAB" will instead be "85,95,171")
#'
#' @return A named vector of colour codes for lymphgen classes and pathology.
#'
#' @import dplyr tidyr
#' @export
#'
#' @examples
#' lymphgen_cols <- get_gambl_colours("lymphgen")
#'
get_gambl_colours <- function(classification = "all",
                              alpha = 1,
                              as_list = FALSE,
                              as_dataframe = FALSE,
                              return_available = FALSE,
                              verbose = FALSE,
                             as_rgb_string=FALSE) {
  all_colours <- list()
  everything <- c()

  #Same as the colours used in IGV
  all_colours[["chromosome"]] =
               c("chr1"="#555FAB",
                "chr2"="#CE3D31",
                "chr3"="#749B58",
                "chr4"="#F0E584",
                "chr5"="#476A85",
                "chr6"="#BA6338",
                "chr7"="#5CB1DD",
                "chr8"="#7F2368",
                "chr9"="#77C269",
                "chr10"="#D595A6",
                "chr11"="#934823",
                "chr12"="#857B8E",
                "chr13"="#C85328",
                "chr14"="#D58F5B",
                "chr15"="#7A65A7",
                "chr16"="#E3AF69",
                "chr17"="#3C1C54",
                "chr18"="#CEDEB7",
                "chr19"="#612B79",
                "chr20"="#AF2064",
                "chr21"="#E6C66F",
                "chr22"="#5B665D",
                "chrX"="#CA992C")

  blood_cols <- c(
    Red = "#c41230",
    Blue = "#115284",
    "Light Green" = "#39b54b",
    Purple = "#5c266c",
    Orange = "#fe9003",
    Green = "#046852",
    Lavendar = "#8781bd",
    "Steel Blue" = "#455564",
    "Light Blue" = "#2cace3",
    Magenta = "#e90c8b",
    LimeGreen = "#a4bb87",
    Brown = "#5f3a17",
    Gray = "#bdbdc1",
    Yellow = "#f9bd1f",
    Mustard = "#b76d29"
  )

  all_colours[["seq_type"]] <- c(
    "mrna" = "#E41A1C",
    "genome" = "#377EB8",
    "capture" = "#4DAF4A"
  )

  all_colours[["type"]] <- c(
    "gain" = "#0000FF",
    "loss" = "#FF0000"
  )

  all_colours[["hmrn"]] <- c(
    "BCL2-MYC" = "#52000F",
    "BCL2" = "#721F0F",
    "SOCS1/SGK1" = "#D66B1F",
    "TET2/SGK1" = "#C41230",
    "MYD88" = "#3B5FAC",
    "NOTCH2" = "#7F3293",
    "NOTCH1" = "#55B55E",
    "Other" = "#ACADAF"
  )

  all_colours[["EBV"]] <- c(
    "EBV-positive" = "#7F055F",
    "EBV-negative" = "#E5A4CB",
    "POS" = "#7F055F",
    "NEG" = "#E5A4CB"
  )

  all_colours[["BL"]] <- c(
    "Q53-BL" = "#A6CEE3",
    "M53-BL" = "#A6CEE3", # added because genetic subgroup still refers to it this way
    "DLBCL-A" = "#721F0F",
    "IC-BL" = "#45425A",
    "DGG-BL" = "#E90C8B",
    "DLBCL-B" = "#FB9A99",
    "DLBCL-C" = "#C41230",
    "DLBCLesque" = "#721F0F",
    "DLBCL-like" = "#721F0F")

  all_colours[["FL"]] <- c(dFL = "#99C1B9", cFL = "#D16666", DLBCL = "#479450",
                          "MEM-like"="#FFB61F", "GC-like"= "#008AEC")

  all_colours[["lymphgenerator"]] <- c(
    "MP3" = "#5B8565",
    "EGB" = "#98622A",
    "ETB" = "#813F3D",
    "aSCI" = "#D66B1F",
    "aSEL" = "#6A0D18",
    "MCaP" = "#5F8CFF",
    "BNZ" = "#8870B6",
    "EZB" = "#721F0F",
    "ST2" = "#C41230",
    "UNCLASS" = "#05631E"
  )

  all_colours[["seq_type"]] <- c(
    "mrna" = "#E41A1C",
    "genome" = "#377EB8",
    "capture" = "#4DAF4A"
  )

  all_colours[["type"]] <- c(
    "gain" = "#0000FF",
    "loss" = "#FF0000"
  )

  all_colours[["hmrn"]] <- c(
    "BCL2-MYC" = "#52000F",
    "BCL2" = "#721F0F",
    "SOCS1/SGK1" = "#D66B1F",
    "TET2/SGK1" = "#C41230",
    "MYD88" = "#3B5FAC",
    "NOTCH2" = "#7F3293",
    "NOTCH1" = "#55B55E",
    "Other" = "#ACADAF"
  )

  all_colours[["EBV"]] <- c(
    "EBV-positive" = "#7F055F",
    "EBV-negative" = "#E5A4CB",
    "POS" = "#7F055F",
    "NEG" = "#E5A4CB"
  )

  all_colours[["BL"]] <- c(
    "Q53-BL" = "#A6CEE3",
    "M53-BL" = "#A6CEE3", # added because genetic subgroup still refers to it this way
    "DLBCL-A" = "#721F0F",
    "IC-BL" = "#45425A",
    "DGG-BL" = "#E90C8B",
    "DLBCL-B" = "#FB9A99",
    "DLBCL-C" = "#C41230"
  )

  all_colours[["FL"]] <- c(dFL = "#99C1B9", cFL = "#D16666", DLBCL = "#479450")

  all_colours[["lymphgenerator"]] <- c(
    "MP3" = "#5B8565",
    "EGB" = "#98622A",
    "ETB" = "#813F3D",
    "aSCI" = "#D66B1F",
    "aSEL" = "#6A0D18",
    "MCaP" = "#5F8CFF",
    "BNZ" = "#8870B6",
    "EZB" = "#721F0F",
    "ST2" = "#C41230",
    "UNCLASS" = "#05631E"
  )

  all_colours[["chapuy_classifier"]] <- c(
    C0 = "#bebebe",
    C1 = "#803D99",
    C2 = "#00A2D2",
    C3 = "#F39123",
    C4 = "#50BFAD",
    C5 = "#DE292A"
  )

  all_colours[["lacy_classifier"]] <- all_colours[["hmrn"]]

  all_colours[["lymphgen"]] <- c(
    "EZB-MYC" = "#52000F",
    "EZB" = "#721F0F",
    "EZB-COMP" = "#C7371A",
    "ST2" = "#C41230",
    "ST2-COMP" = "#EC3251",
    "MCD" = "#3B5FAC",
    "MCD-COMP" = "#6787CB",
    "BN2" = "#7F3293",
    "BN2-COMP" = "#A949C1",
    "N1" = "#55B55E",
    "N1-COMP" = "#7FC787",
    "A53" = "#5b6d8a",
    "Other" = "#ACADAF",
    "COMPOSITE" = "#ACADAF"
  )


  all_colours[["mutation"]] <-
    c(
      "Nonsense_Mutation" = unname(blood_cols["Red"]),
      "Missense_Mutation" = unname(blood_cols["Light Green"]),
      "Multi_Hit" = unname(blood_cols["Steel Blue"]),
      "Frame_Shift_Ins" = unname(blood_cols["Magenta"]),
      "Frame_Shift_Del" = unname(blood_cols["Magenta"]),
      "In_Frame_Ins" = unname(blood_cols["Brown"]),
      "In_Frame_Del" = unname(blood_cols["Brown"]),
      "Nonstop_Mutation" = unname(blood_cols["Light Blue"]),
      "Translation_Start_Site" = unname(blood_cols["Lavendar"]),
      "Splice_Site" = unname(blood_cols["Orange"]),
      "Splice_Region" = unname(blood_cols["Orange"]),
      "3'UTR" = unname(blood_cols["Yellow"]),
      "5'UTR" = unname(blood_cols["LimeGreen"]),
      "Intron" = unname(blood_cols["Mustard"]),
      "Silent" = "#D8A7CA"
    )

  all_colours[["rainfall"]] <-
    c(
      "C>A" = "#2196F3FF",
      "C>G" = "#3F51B5FF",
      "C>T" = "#F44336FF",
      "InDel" = "#A020F0",
      "T>A" = "#4CAF50FF",
      "T>C" = "#FFC107FF",
      "T>G" = "#FF9800FF"
    )

  all_colours[["pos_neg"]] <- c(
    "POS" = "#c41230",
    "NEG" = "#E88873",
    "PARTIAL" = "#E88873",
    "yes" = "#c41230",
    "no" = "#E88873",
    "YES" = "#c41230",
    "NO" = "#E88873",
    "FAIL" = "#bdbdc1",
    "positive" = "#c41230",
    "negative" = "#E88873",
    "1" = "#c41230",
    "0" = "#E8887366",
    "fail" = "#bdbdc1"
  )

  all_colours[["copy_number"]] <- c(
    "nLOH" = "#E026D7",
    "14" = "#380015",
    "15" = "#380015",
    "13" = "#380015",
    "12" = "#380015",
    "11" = "#380015",
    "10" = "#380015",
    "9" = "#380015",
    "8" = "#380015",
    "7" = "#380015",
    "6" = "#380015",
    "5" = "#67001F",
    "4" = "#B2182B",
    "3" = "#D6604D",
    "2" = "#ede4c7",
    "1" = "#92C5DE",
    "0" = "#4393C3"
  )

  all_colours[["aneuploidy"]]=c(
    "iso-pq_lossgain"="#E208D7",
    "iso-qp_lossgain"="#E563D7",
    "arm-p_gain"="#6B641F",
    "chrom_gain"="#6C331F",
    "arm-q_gain"="#6A504D",
    "chrom_loss"="#4693C3",
    "arm-p_loss"="#93B6DE",
    "arm-q_loss"="#90E0DE"
  )

  all_colours[["blood"]] <- c(
    "Red" = "#c41230", "Blue" = "#115284", "Light Green" = "#39b54b",
    "Purple" = "#5c266c", "Orange" = "#fe9003", "Green" = "#046852",
    "Lavendar" = "#8781bd", "Steel Blue" = "#455564",
    "Light Blue" = "#2cace3", "Magenta" = "#e90c8b", "Mustard" = "#b76d29",
    "LimeGreen" = "#a4bb87", "Brown" = "#5f3a17", "Gray" = "#bdbdc1",
    "Yellow" = "#f9bd1f"
  )
  all_colours[["sex"]] <- c(
    "M" = "#118AB2",
    "Male" = "#118AB2",
    "male" = "#118AB2",
    "F" = "#EF476F",
    "Female" = "#EF476F",
    "female" = "#EF476F"
  )

  all_colours[["clinical"]] <-
    c(
      "M" = "#118AB2",
      "Male" = "#118AB2",
      "F" = "#EF476F",
      "Female" = "#EF476F",
      "EBV-positive" = "#7F055F",
      "EBV-negative" = "#E5A4CB",
      "POS" = "#c41230",
      "NEG" = "#E88873",
      "FAIL" = "#bdbdc1",
      "Alive" = "#046852",
      "alive" = "#046852",
      "dead" = "#a4bb87",
      "Dead" = "#a4bb87",
      "deceased" = "#a4bb87",
      "unknown" = "#C3C9E9",
      "IPI_0" = "#3B9AB2",
      "IPI_1" = "#78B7C5",
      "IPI_2" = "#EBCC2A",
      "IPI_3" = "#E1AF00",
      "IPI_4" = "#F21A00",
      "Adult" = "#DCE0E5",
      "adult" = "#DCE0E5",
      "Pediatric" = "#677A8E",
      "pediatric" = "#677A8E",
      "Diagnosis" = "#E57A44",
      "A" = "#E57A44",
      "B" = "#721817",
      "C" = "#721817",
      "D" = "#721817",
      "E" = "#721817",
      "Progression" = "#A44A3F",
      "Relapse" = "#721817",
      "I" = "#75F4F4",
      "FOLL1" = "#75F4F4",
      "II" = "#90E0F3",
      "FOLL2" = "#90E0F3",
      "IIIA" = "#B8B3E9",
      "FOLL3A" = "#B8B3E9",
      "IIIB" = "#D999B9",
      "FOLL3B" = "#D999B9",
      "matched" = "#F0B67F",
      "unmatched" = "#D6D1B1",
      "FF" = "#009FFD",
      "frozen" = "#009FFD",
      "FFPE" = "#95B2B8",
      "ctDNA" = "#7E6148",
      "NA" = "white"
    )
  all_colours[["pathology"]] <- c(
    "B-ALL" = "#C1C64B",
    "CLL" = "#889BE5",
    "MCL" = "#40E0D0",
    "BL" = "#926CAD",
    "mBL" = "#34C7F4",
    "tFL" = "#FF8595",
    "DLBCL-BL-like" = "#34C7F4",
    "pre-HT" = "#754F5B",
    "PMBL" = "#227C9D",
    "PMBCL" = "#227C9D",
    "FL" = "#EA8368",
    "no-HT" = "#EA8368",
    "COMFL" = "#8BBC98",
    "COM" = "#8BBC98",
    "post-HT" = "#479450",
    "DLBCL" = "#479450",
    "denovo-DLBCL" = "#479450",
    "HGBL-NOS" = "#294936",
    "HGBL" = "#294936",
    "HGBL-DH/TH" = "#7A1616",
    "PBL" = "#E058C0",
    "Plasmablastic" = "#E058C0",
    "CNS" = "#E2EF60",
    "cHL"="#C1C15C",
    "THRLBCL" = "#A5F2B3",
    "MM" = "#CC9A42",
    "SCBC" = "#8c9c90",
    "UNSPECIFIED" = "#cfba7c",
    "OTHER" = "#cfba7c",
    "MZL" = "#065A7F",
    "SMZL" = "#065A7F",
    "Prolymphocytic" = "#7842f5"
  )
  all_colours[["coo"]] <- c(
    "ABC" = "#05ACEF",
    "UNCLASS" = "#05631E",
    "Unclass" = "#05631E",
    "U" = "#05631E",
    "UNC" = "#05631E",
    "GCB" = "#F58F20",
    "DHITsig-" = "#F58F20",
    "DHITsigNeg" = "#F58F20",
    "DHITsig-IND" = "#003049",
    "DHITsig+" = "#D62828",
    "DHITsigPos" = "#D62828",
    "NA" = "#ACADAF"
  )
  all_colours[["cohort"]] <- c(
    "Chapuy" = "#8B0000", "Chapuy, 2018" = "#8B0000",
    "Arthur" = "#8845A8", "Arthur, 2018" = "#8845A8",
    "Schmitz" = "#2C72B2", "Schmitz, 2018" = "#2C72B2",
    "Reddy" = "#E561C3", "Reddy, 2017" = "#E561C3",
    "Morin" = "#8DB753", "Morin, 2013" = "#8DB753",
    "Kridel" = "#4686B7", "Kridel, 2016" = "#4686B7",
    "ICGC" = "#E09C3B", "ICGC, 2018" = "#E09C3B",
    "Grande" = "#e90c8b", "Grande, 2019" = "#e90c8b"
  )

  all_colours[["indels"]] <- c("DEL" = "#53B1FC", "INS" = "#FC9C6D")
  all_colours[["svs"]] <- c("DEL" = "#53B1FC", "DUP" = "#FC9C6D")
  all_colours[["genetic_subgroup"]] <- c(all_colours[["lymphgen"]], all_colours[["BL"]], all_colours[["FL"]])
  all_colours[["domains"]] <- c(
    "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
    "#084594", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A",
    "#E31A1C", "#B10026", "#91003F", "#CE1256", "#E7298A", "#DF65B0", "#C994C7",
    "#D4B9DA", "#E7E1EF", "#F7F4F9"
  )
  names(all_colours[["domains"]]) <- paste0(
    "domain",
    1:length(all_colours[["domains"]])
  )

  # print(all_colours)
  if (alpha < 1) {
    for (colslot in names(all_colours)) {
      raw_cols <- all_colours[[colslot]]
      raw_cols_rgb <- col2rgb(raw_cols)
      alpha_cols <- rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), maxColorValue = 255L)
      names(alpha_cols) <- names(raw_cols)
      all_colours[[colslot]] <- alpha_cols
    }
  }

  if(as_rgb_string){
    for(colslot in names(all_colours)){
      raw_cols = all_colours[[colslot]]
      raw_cols_rgb = col2rgb(raw_cols)
      col_string = paste(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],sep=",")
      names(col_string) = names(raw_cols)
      all_colours[[colslot]] = col_string
    }
  }

  for (this_group in names(all_colours)) {
    everything <- c(everything, all_colours[[this_group]])
  }
  # return matching value from lowercase version of the argument if it exists
  lc_class <- tolower(classification)
  if (return_available) {
    return(names(all_colours))
  }
  if (classification %in% names(all_colours)) {
    if (as_dataframe) {
      some_col <- all_colours[[classification]]
      df_ugly <- data.frame(name = names(some_col), colour = unname(some_col))
      df_tidy <- mutate(df_ugly, group = classification)
      return(df_tidy)
    }
    return(all_colours[[classification]])
  } else if (lc_class %in% names(all_colours)) {
    return(all_colours[[lc_class]])
  } else if (as_list) {
    return(all_colours)
  } else if (as_dataframe) {
    df_ugly <- data.frame(name = names(unlist(all_colours, use.names = T)), colour = unlist(all_colours, use.names = T))
    df_tidy <- separate(df_ugly, name, into = c("group", "name"), sep = "\\.")
    return(df_tidy)
  } else {
    return(everything)
  }
}
