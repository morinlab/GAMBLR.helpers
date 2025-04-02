#' @title Trim Scale Expressions.
#'
#' @description INTERNAL HELPER FUNCTION called by prettyOncoplot, not meant for
#' out-of-package usage.
#'
#' @details INTERNAL FUNCTION called by prettyOncoplot, not meant for
#' out-of-package usage.
#'
#' @param x Numeric value (of expression) to be trimmed.
#'
#' @return Numeric value.
#'
#' @examples
#' trimmed = trim_scale_expression(c(2,4,6,4,4.5))
#' trimmed
#' 
#' @export
trim_scale_expression = function(x){
  quants = unname(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  x = ifelse(x < quants[1], quants[1], x)
  x = ifelse(x > quants[2], quants[2], x)
  x = (x - quants[1]) / (quants[2] - quants[1])
  return(x)
}
