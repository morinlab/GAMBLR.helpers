#' @title Cool overlap of data frames.
#'
#' @description This function implements overlap of 2 data frames that contain
#' regions of coordinates similar to what data.table::foverlaps does. Unlike
#' foverlaps, this function takes as input data frame class objects, and relies
#' on dplyr solution rather than data.table handling, therefore allowing usage
#' of data frames with virtually unlimited dimensions without crashing. This
#' implementation uses same logic of different types of overlaps as the original
#' foverlaps solution (any, start, end, within, equal). The type any is default
#' and allows for any overlapping solution between 2 regions. The type start
#' only considers regions with exact same start position as overlap; similarly
#' type end considers regions overlapped when the end positions are exact
#' matches. Type within means that regions are overlapped when one is contained
#' in another and neither start nor end positions match. Finally, type equal
#' only considers overlap when both start and end positions match for both
#' regions. For any type, the presence of any additional column not directly
#' specifying regions (for example, Chromosome) will serve similar to a grouping
#' variable.
#' The generated output of this function will contain the overlapping regions
#' and all columns present in the data frame data1, as well as any columns from
#' the data frame wupplied with data2 argument, except for those columns present
#' in data2 that are used for overlap.
#'
#' @param data1 Data frame with data to overlap. Required parameter. The minimal
#'      required columns are those supplied with the argument columns1. Will
#'      dictate the naming of the columns used for overlap in the output.
#' @param data2 Data frame with data to overlap. Required parameter. The minimal
#'      required columns are those supplied with the argument columns2.
#' @param columns1 The list of columns from data frame data1 to be used to find
#'      overlapping regions.
#' @param columns2 The list of columns from data frame data2 to be used to find
#'      overlapping regions.
#' @param type Character specifying the way to find overlaps. Accepted values
#'      are any (used as default), start, end, within, and equal. Please see
#'      function description for more details of different types.
#'
#' @return data frame
#'
#' @examples
#' library(GAMBLR.data)
#' maf <- get_coding_ssm()
#' onco_matrix <- create_onco_matrix(maf_df = maf)
#'
#' @import dplyr
#' @export
#'
cool_overlaps <- function(
    data1,
    data2,
    columns1 = c("Chromosome", "Start_Position", "End_Position"),
    columns2 = c("Chromosome", "Start_Position", "End_Position"),
    type = "any"
){

    # Ensure all columns provided for overlap are present in the data frame
    if(! length(columns1) == length(intersect(columns1, colnames(data1)))){
        stop(
            "Not all of the requested columns for overlap in data1 are present."
        )
    }

    if(! length(columns2) == length(intersect(columns2, colnames(data2)))){
        stop(
            "Not all of the requested columns for overlap in data2 are present."
        )
    }

    # What is the name of the column in columns1 that specifies start and end?
    start1 <- columns1[grepl("start", columns1, ignore.case = TRUE)]
    end1 <- columns1[grepl("end", columns1, ignore.case = TRUE)]

    # What is the name of the column in columns1 that specifies start and end?
    start2 <- columns1[grepl("start", columns2, ignore.case = TRUE)]
    end2 <- columns1[grepl("end", columns2, ignore.case = TRUE)]

    # Prepare for overlap
    overlap <- dplyr::inner_join(
        data1,
        data2,
        by = structure(names = columns1, .Data = columns2),
        relationship = "many-to-many"
    )

    # Return matches based on mode
    if(type == "any"){
        message(
            "Running in default mode of any..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(start2) <= !!sym(end1)
            )
    } else if (type == "start"){
        message(
            "Running in the mode start..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(start1) == !!sym(start2)
            )
    } else if (type == "end"){
        message(
            "Running in the mode end..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               !!sym(end1) == !!sym(end2)
            )
    } else if (type == "within"){
        message(
            "Running in the mode within..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               (!!sym(start1) >= !!sym(start2)) & (!!sym(end1) <= !!sym(end2))
            )
    } else if (type == "equal"){
        message(
            "Running in the mode equal..."
        )
        overlap <- overlap %>%
            dplyr::filter(
               (!!sym(start1) == !!sym(start2)) & (!!sym(end1) == !!sym(end2))
            )
    } else {
        message(
            "You have requested mode that is not supported."
        )
        stop(
            "Please supply one of any, start, end, within, or equal with type."
        )
    }

    return(overlap)
}
