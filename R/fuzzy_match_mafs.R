#' @title Fuzzy match mafs
#'
#' @description Compare two MAFs and identify pairs of rows as likely matches
#' using approximate positions (mostly for indels)
#'
#' @details Specify the two mafs in a data frame format to compare as `maf1` and
#' `maf2`. In the generated output, the data frame is derived from maf1 with two
#' new columns (maf1_ID and maf2_ID). Values in maf1_ID are derived from the
#' maf1. The values in maf2_ID are derived from maf2 and are NA for any rows
#' that could not be matched to maf2.
#'
#' @param maf1 The index maf (akin to the first argument in a left_join). The
#'      expected input format is data frame. The minimal required columns are
#'      Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position,
#'      Reference_Allele, and Tumor_Seq_Allele2.
#' @param maf2 The maf you will search for potential matches to rows in maf1
#'      (akin to the second argument in a left_join). Has the same minimal
#'      required columns as the maf file supplied with maf1 argument.
#' @param padding Number of nucleotides of padding to add around Start_Position
#'      and End_Position for fuzzy match. The default is 0 (no padding).
#' @param matching_only Set to TRUE if you only want the rows of maf1 that were
#'      successfully matched to maf2. Default is FALSE.
#' @param invert Set to TRUE if you ony want the rows of maf1 that could not be
#'      matched to maf2. Default is FALSE.
#'
#' @return data frame
#'
#' @import dplyr tidyr tibble
#' @export
#'
#' @examples
#' maf_1 <- get_ssm_by_sample(these_sample_id = "DOHH-2") %>%
#'     dplyr::filter(Chromosome == "1")
#'
#' maf_2 <- get_ssm_by_sample(these_sample_id = "SU-DHL-4") %>%
#'     dplyr::filter(Chromosome == "1")
#'
#' matched_maf <- fuzzy_match_mafs(
#'     maf1 = maf_1,
#'     maf2 = maf_2
#' )
#'
fuzzy_match_mafs <- function(
    maf1,
    maf2,
    padding = 0,
    matching_only = FALSE,
    invert = FALSE
){
    maf1 <- maf1 %>%
        as.data.frame()
    maf2 <- maf2 %>%
        dplyr::select(
            Tumor_Sample_Barcode,
            Chromosome,
            Start_Position,
            End_Position,
            Reference_Allele,
            Tumor_Seq_Allele2
        ) %>%
        dplyr::mutate(
            Start_Position = Start_Position - padding,
            End_Position = End_Position + padding
        ) %>%
        as.data.frame()

    maf1 <- dplyr::mutate(
        maf1,
        maf1_ID = paste(
            Tumor_Sample_Barcode,
            Chromosome,
            Start_Position,
            End_Position,
            Reference_Allele,
            Tumor_Seq_Allele2,
            sep = "_"
        )
    )

    maf2 <- dplyr::mutate(
        maf2,
        maf2_ID = paste(
            Tumor_Sample_Barcode,
            Chromosome,
            Start_Position,
            End_Position,
            Reference_Allele,
            Tumor_Seq_Allele2,
            sep = "_"
        )
    )

    # naive matching: overlap within Start_Position - End_Position
    columns_to_overlap = c(
        "Chromosome",
        "Tumor_Sample_Barcode",
        "Start_Position",
        "End_Position"
    )
    matched <- cool_overlaps(
        data1 = maf1,
        data2 = maf2,
        columns1 = columns_to_overlap,
        columns2 = columns_to_overlap,
        type = "any"
    ) %>%
        unique() %>%
        as.data.frame()

    if(matching_only) {
        matched <- dplyr::filter(matched, !is.na(maf2_ID))
    } else if(invert) {
        matched <- dplyr::filter(matched, is.na(maf2_ID))
    }

    matched <- dplyr::select(
        matched,
        -contains(".y")
    ) %>%
    dplyr::rename(
        c(
            "Start_Position" = "Start_Position.x",
            "Tumor_Seq_Allele2" = "Tumor_Seq_Allele2.x",
            "Reference_Allele" = "Reference_Allele.x",
            "End_Position" = "End_Position.x"
        )
    )

    return(matched)
}
