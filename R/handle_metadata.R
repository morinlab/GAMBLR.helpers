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
#' @import GAMBLR.data
#'
#' @noRd
#'
#' @export
handle_metadata = function() {
    if ("GAMBLR.data" %in% installed.packages()) {
        return(GAMBLR.data::gambl_metadata)
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
