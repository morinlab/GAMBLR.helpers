#' @title Create onco matrix from maf data.
#'
#' @description This function implements generation of ComplexHeatmap-compatible
#' matrix using dplyr. The output of this function is identical to the implementation
#' in maftools, hovewer the main advantage of this version is tidy syntax and
#' avoidance of data.table dependency since the latter package cannot handle
#' large tables required for processing of large multi-sample maf files.
#'
#' @param maf_df Data frame with maf data. Required parameter. The minimal required columns are Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Start_Position, and End_Position.
#' @param genes List of genes to return in the resulting matrix. When not provided, matrix is generated for each gene present in the input maf data.
#' @param add_missing When list of genes is provided and some of the specified genes are not mutated in the provided maf data, they will still be added to the matrix with 0% mutation frequency. Default TRUE.
#'
#' @return matrix
#'
#' @examples
#' library(GAMBLR.data)
#' maf <- get_coding_ssm()
#' onco_matrix <- create_onco_matrix(maf_df = maf)
#'
#' @import dplyr tidyr tibble
#' @export
#'
create_onco_matrix = function(
    maf_df,
    genes,
    add_missing = TRUE
){
    if(missing(maf_df)){
        stop(
            "Data frame with maf data was not provided."
        )
    }

    if(! missing(genes)){
        maf_df <- maf_df %>%
            dplyr::filter(
                Hugo_Symbol %in% genes
            )
    }

    onco_matrix_coding <- coding_class[
        !coding_class %in% c("Silent", "Splice_Region", "Targeted_Region")
    ]


    onco_matrix <- maf_df %>%
        dplyr::distinct(
            Tumor_Sample_Barcode, Hugo_Symbol,
            Variant_Classification,
            Start_Position, End_Position
        ) %>%
        dplyr::select(
            Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification
        ) %>%
        dplyr::filter(
            Variant_Classification %in% onco_matrix_coding
        ) %>%
        dplyr::group_by(
            Hugo_Symbol, Tumor_Sample_Barcode
        ) %>%
        dplyr::mutate(
            n = n(),
            Variant_Classification = ifelse(
                n > 1,
                "Multi_Hit",
                Variant_Classification
            )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::distinct(
            Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification
        ) %>%
        tidyr::pivot_wider(
            names_from = "Tumor_Sample_Barcode",
            values_from = "Variant_Classification"
        ) %>%
        tibble::column_to_rownames("Hugo_Symbol") %>%
        replace(is.na(.), "") %>%
        as.matrix

    if(!missing(genes) & add_missing){
        missing_genes <- setdiff(
            genes,
            rownames(onco_matrix)
        )
        if (length(missing_genes) == 0){
            message("There are no missing genes during onco matrix assembly.")
        } else {
            message(
                paste0(
                    "Adding ",
                    length(missing_genes),
                    " gene(s) with 0% mut frequency during onco matrix assembly."
                )
            )
            missing_matrix <- matrix(
                    data = rep("", ncol(onco_matrix)*length(missing_genes)),
                    nrow = length(missing_genes),
                    ncol = ncol(onco_matrix),
                    byrow = FALSE,
                    dimnames = list(missing_genes, colnames(onco_matrix))
                )

            onco_matrix <- rbind(
                onco_matrix,
                missing_matrix
            )
        }
    }


    return(onco_matrix)
}
