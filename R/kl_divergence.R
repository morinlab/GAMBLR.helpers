#' Calculate Kullback-Leibler Divergence
#'
#' This function computes the Kullback-Leibler (KL) divergence between two
#'      probability distributions, with an optional small constant (epsilon)
#'      added to avoid zero probabilities, which would otherwise cause division
#'      by zero or undefined logarithms.
#'
#' @param P A numeric vector representing the first probability distribution.
#'      The sum of "P" should be 1, but the function will normalize it if
#'      necessary.
#' @param Q A numeric vector representing the second probability distribution.
#'      The sum of "Q" should be 1, but the function will normalize it if
#'      necessary.
#' @param epsilon A small positive number (default = 1e-7) to be added to each
#'      probability in P and Q to avoid zero probabilities. This helps to
#'      prevent division by zero or log(0).
#'
#' @return float
#'
#' @examples
#' P <- c(0.1, 0.4, 0.3, 0.2)
#' Q <- c(0.2, 0.3, 0.4, 0.1)
#'
#' kl_divergence(P, Q)
#'
#' @export
kl_divergence <- function(P, Q, epsilon = 1e-7) {
    P <- P + epsilon
    Q <- Q + epsilon

    P <- P / sum(P)
    Q <- Q / sum(Q)

    # KL divergence formula: sum(P * log(P / Q))
    return(sum(P * log(P / Q), na.rm = TRUE))
}
