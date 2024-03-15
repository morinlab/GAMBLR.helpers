#' @title Copy not clobber.
#'
#' @description TODO.
#'
#' @param from_file TODO.
#' @param to_file TODO.
#' @param force TODO.
#'
#' @return full path to the file that was written
#'
#' @export
copy_no_clobber = function(from_file,
                           to_file,
                           force = FALSE){

  to_dir = dirname(to_file)
  suppressMessages(suppressWarnings(dir.create(to_dir,recursive = T)))
  print(paste("COPYING",from_file,"TO",to_file))
  if(force){
    file.copy(from_file,to_file)
  }
  return(to_file)
}
