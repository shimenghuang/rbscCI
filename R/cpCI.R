#' Clopper-Pearson confidence interval
#'
#' @param n_tot Total number of experiments
#' @param n_suc Number of successes
#' @param conf Confidence level (1-alpha)
#' @details ...
#'
#' @export cpCI
cpCI <- function(n_tot, n_suc, conf) {
  .cpCI(n_tot, n_suc, conf)
}
