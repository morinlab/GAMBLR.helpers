#' @title Handle metadata from bundled packages.
#'
#' @description Will return the metadata from the GAMBLR.data package when it is
#' installed. Otherwise, will instruct user to call
#' GAMBLR.results::get_gambl_metadata().
#'
#' @details INTERNAL FUNCTION for handling the metadata.
#'
#' @return data frame
#'
#' @param this_seq_type The seq type of the samples to return with metadata.
#'
#' @import GAMBLR.data dplyr
#'
#' @export
handle_metadata = function(this_seq_type = "genome") {
    if ("GAMBLR.data" %in% installed.packages()) {
        message(
            "Using bundled metadata from GAMBLR.data"
        )
        return(
            GAMBLR.data::gambl_metadata %>%
                dplyr::filter(seq_type %in% this_seq_type)
        )
    }else if ("GAMBLR.results" %in% installed.packages()) {
        stop(
            "Use GAMBLR.results::get_gambl_metadata() to obtain the metadata."
        )
    }else {
        stop(
            "Please install GAMBLR.data or GAMBLR.results to obtain metadata."
        )
    }
}
