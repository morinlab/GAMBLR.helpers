#' Check and Clean Metadata
#'
#' This function validates and optionally cleans a metadata data frame so that it
#' conforms to the expected format for downstream functions. Specifically, it checks that:
#'
#' \itemize{
#'   \item A \code{sample_id} column exists and contains no \code{NA} values.
#'   \item A \code{seq_type} column exists and all its values are one of \code{"mrna"},
#'         \code{"genome"}, or \code{"capture"}.
#'   \item A \code{pathology} column exists.
#'   \item No two rows have identical combinations of \code{sample_id} and \code{seq_type}.
#' }
#'
#' The user can specify how to handle each problem:
#'
#' \itemize{
#'   \item \strong{NA in \code{sample_id}:} Set \code{na_sample_id_action} to \code{"stop"}
#'         (default) to throw an error or to \code{"drop"} to remove offending rows.
#'   \item \strong{Invalid \code{seq_type}:} Set \code{invalid_seq_type_action} to \code{"stop"}
#'         (default) to throw an error or to \code{"drop"} to remove rows with invalid values.
#'   \item \strong{Duplicate combinations of \code{sample_id} and \code{seq_type}:}
#'         Set \code{duplicate_action} to \code{"stop"} (default), \code{"keep_first"}, or \code{"keep_last"}.
#' }
#'
#' @param df A data frame containing the metadata.
#' @param na_sample_id_action How to handle \code{NA} values in \code{sample_id}.
#'        Options are \code{"stop"} (default) or \code{"drop"}.
#' @param invalid_seq_type_action How to handle invalid values in \code{seq_type}.
#'        Options are \code{"stop"} (default) or \code{"drop"}.
#' @param duplicate_action How to handle duplicate combinations of \code{sample_id} and \code{seq_type}.
#'        Options are \code{"stop"} (default), \code{"keep_first"}, or \code{"keep_last"}.
#'
#' @return A cleaned data frame that meets the expected metadata requirements.
#' @import dplyr
#' @export
check_and_clean_metadata <- function(df,
                                     na_sample_id_action = c("stop", "drop"),
                                     invalid_seq_type_action = c("stop", "drop"),
                                     duplicate_action = c("stop", "keep_first", "keep_last")) {
  # Match the argument options
  na_sample_id_action <- match.arg(na_sample_id_action)
  invalid_seq_type_action <- match.arg(invalid_seq_type_action)
  duplicate_action <- match.arg(duplicate_action)
  
  # Validate required columns exist using base R
  required_cols <- c("sample_id", "seq_type", "pathology")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing: ", paste(missing_cols, collapse = ", "))
  }
  
  # Use dplyr to remove rows with NA in sample_id if requested
  if (any(is.na(df$sample_id))) {
    if (na_sample_id_action == "stop") {
      stop("The 'sample_id' column contains NA values.")
    } else if (na_sample_id_action == "drop") {
      df <- df %>% filter(!is.na(sample_id))
      message("Rows with NA in 'sample_id' have been dropped.")
    }
  }
  
  # Define allowed values for seq_type
  allowed_seq <- c("mrna", "genome", "capture")
  invalid_seq <- df %>% 
    filter(!seq_type %in% allowed_seq) %>% 
    distinct(seq_type) %>% 
    pull(seq_type)
  if (length(invalid_seq) > 0) {
    if (invalid_seq_type_action == "stop") {
      stop("The 'seq_type' column contains invalid values: ", paste(invalid_seq, collapse = ", "))
    } else if (invalid_seq_type_action == "drop") {
      df <- df %>% filter(seq_type %in% allowed_seq)
      message("Rows with invalid 'seq_type' values (", paste(invalid_seq, collapse = ", "), ") have been dropped.")
    }
  }
  
  # Resolve duplicate combinations of sample_id and seq_type
  dup_summary <- df %>%
    group_by(sample_id, seq_type) %>%
    tally() %>%
    filter(n > 1)
  
  if (nrow(dup_summary) > 0) {
    if (duplicate_action == "stop") {
      stop("Duplicate combinations of 'sample_id' and 'seq_type' were found.")
    } else if (duplicate_action == "keep_first") {
      df <- df %>% 
        group_by(sample_id, seq_type) %>% 
        slice(1) %>% 
        ungroup()
      message("Duplicate rows (keeping first occurrence) for 'sample_id' and 'seq_type' have been dropped.")
    } else if (duplicate_action == "keep_last") {
      df <- df %>% 
        group_by(sample_id, seq_type) %>% 
        slice_tail(n = 1) %>% 
        ungroup()
      message("Duplicate rows (keeping last occurrence) for 'sample_id' and 'seq_type' have been dropped.")
    }
  }
  
  return(df)
}
