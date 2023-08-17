#' @export
socketWrite = function(sock, string){
  print(string)
  write.socket(sock, string)
  response <- read.socket(sock)
  return(response)
}
