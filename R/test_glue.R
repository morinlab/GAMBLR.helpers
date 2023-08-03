#' @export
test_glue = function(placeholder="INSERTED"){
  some_string = "this text has {placeholder}"
  print(some_string)
  ss=glue::glue(some_string)
  print(ss)
}
