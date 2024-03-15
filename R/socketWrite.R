#' @title Socket write.
#'
#' @description TODO.
#'
#' @param sock TODO.
#' @param string TODO.
#'
#' @return someting
#'
#' @export
socketWrite = function(sock, string){
  print(string)
  write.socket(sock, string)
  response <- read.socket(sock)
  return(response)
}
