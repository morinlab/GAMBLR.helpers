% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_expression_data.R
\name{normalize_expression_data}
\alias{normalize_expression_data}
\title{Normalize Expression Data}
\usage{
normalize_expression_data(
  expression_df,
  log_transform = FALSE,
  exclude_zero_negative_values = FALSE
)
}
\arguments{
\item{expression_df}{Input data frame with raw gene expression counts. This is a required parameter.}

\item{log_transform}{Logical. Should the data be log-transformed before normalization? Default is FALSE.}

\item{exclude_zero_negative_values}{Logical. Should zero/negative values be excluded when normalizing the data? Default is FALSE.}
}
\value{
A data frame containing the normalized expression data with z-scores.
}
\description{
This function normalizes raw gene expression data, calculates cBioPortal
styled z-scores, and optionally log-transforms the data and excludes zero and
negative values. This implementation uses the whole population
of samples in the provided data to perform normalization. Each gene is
normalized separately. The expression distribution of the gene is estimated
by calculating the mean and variance of the expression values for all samples
in the reference poplulation.
}
\details{
The input data frame is specified as expression_df parameter. The expression
data must contain gene identifier as first column, followed by the data
columns for each individual sample. This function is agnostic to the
identifier column name (hugo symbol or ENSG ids) and expects the data to
contain more than one unique sample in order to properly calculate the
z-score. Each row in the input data frame represents feature (gene). The
values should represent read counts or RPKM/FPKM for the RNA-Seq data. Please
refer to the cBioPortal documentation for more information on the gene
expression data requirements.
}
\examples{
\dontrun{
# Normalize expression data with log transformation and exclusion of zero/negative values
normalize_expression_data(expression_df, log_transform = TRUE, exclude_zero_negative_values = TRUE)
}

}
