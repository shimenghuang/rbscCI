#' Blyth-Still-Casella confidence interval
#'
#' @param n_tot Total number of experiments
#' @param n_suc Number of successes
#' @param conf Confidence level (1-alpha)
#' @details ...
#'
#' @export bscCI
bscCI <- function(n_tot, n_suc, conf) {
  .bscCI(n_tot, n_suc, conf)
}
