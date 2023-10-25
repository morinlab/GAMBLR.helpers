#' Normalize Expression Data
#'
#' This function normalizes raw gene expression data, calculates cBioPortal
#' styled z-scores, and optionally log-transforms the data and excludes zero and
#' negative values. This implementation uses the whole population
#' of samples in the provided data to perform normalization. Each gene is
#' normalized separately. The expression distribution of the gene is estimated
#' by calculating the mean and variance of the expression values for all samples
#' in the reference poplulation.
#' 
#' The input data frame is specified as expression_df parameter. The expression
#' data must contain gene identifier as first column, followed by the data
#' columns for each individual sample. This function is agnostic to the
#' identifier column name (hugo symbol or ENSG ids) and expects the data to
#' contain more than one unique sample in order to properly calculate the
#' z-score. Each row in the input data frame represents feature (gene). The
#' values should represent read counts or RPKM/FPKM for the RNA-Seq data. Please 
#' refer to the cBioPortal documentation for more information on the gene
#' expression data requirements.
#'
#' @param expression_df Input data frame with raw gene expression counts. This is a required parameter.
#' @param log_transform Logical. Should the data be log-transformed before normalization? Default is FALSE.
#' @param exclude_zero_negative_values Logical. Should zero/negative values be excluded when normalizing the data? Default is FALSE.
#'
#' @return A data frame containing the normalized expression data with z-scores.
#' @export
#'
#' @examples
#' \dontrun{
#' # Normalize expression data with log transformation and exclusion of zero/negative values
#' normalize_expression_data(expression_df, log_transform = TRUE, exclude_zero_negative_values = TRUE)
#' }
#'
#' @import dplyr readr
#'

normalize_expression_data <- function(
    expression_df,
    log_transform = FALSE,
    exclude_zero_negative_values = FALSE
) {

    if(missing(expression_df)){
        stop(
            "No expression data frame is provided."
        )
    }

    if (exclude_zero_negative_values) {
        message(
            "Zeros and negative counts will be excluded from the reference population when calculating the zscores."
        )
    }

    # How many samples are in the data
    sample_n <- length(
        unique(
            colnames(
                expression_df %>%
                    select(where(is.numeric))
            )
        )
    )
    # Make sure there is more than one sample
    if (sample_n <= 1) {
        message("Cannot calculate z-scores")
        stop("Expression data contains one or no samples.")
    }


    # Extract and transform expression data
    expression_data <- expression_df %>%
        select(where(is.numeric))

    if (log_transform) {
        message(
            "Performing log-transformation..."
        )
        expression_data <- log2(expression_data + 1)
    }


    # Some helpers

    # Calculate mean and std for the samples whose values are not Null or NA (n)
    # for rnaseq data, ignore the negative and zero values.
    calculate_mean_std <- function(data, exclude_zero_negative_values) {
        # Drop NA
        filtered_data <- data[!is.na(data)]
        # Drop 0 and negative
        if (exclude_zero_negative_values){
            filtered_data <- data[data > 0]
        }

        n <- length(filtered_data)

        if (n <= 1) {
            mu <- 0
            sigma <- 0
        } else {
            mu <- mean(filtered_data, na.rm = TRUE)
            sigma <- sd(filtered_data, na.rm = TRUE)
        }
        return(c(mu, sigma))
    }

    # Calculate z-scores
    mu_sigma <- expression_data %>%
        rowwise() %>%
        mutate(
            mu = calculate_mean_std(
                c_across(where(is.numeric)),
                exclude_zero_negative_values
            )[1],
            sigma = calculate_mean_std(
                c_across(where(is.numeric)),
                exclude_zero_negative_values
            )[2]
        )

    zscores <- expression_df %>%
        select(!where(is.numeric)) %>%
        bind_cols(
            mu_sigma %>%
                ungroup() %>% 
                mutate(across(1:(ncol(.)-2), ~ (. - mu)/sigma)) %>%
                select( - c(mu, sigma))
        )

    return(zscores)

}
