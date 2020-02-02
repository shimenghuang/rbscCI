#' Blyth-Still-Casella confidence interval
#' 
#' @details Comuputes the exact Blyth-Still-Casella binomial confidence interval. The initial CI is the Clopper-Pearson confidence interval.
#'
#' @param n_tot Total number of experiments
#' @param n_suc Number of successes
#' @param conf Confidence level (1-alpha)
#' @details ...
#'
#' @export bscCI
bscCI <- function(n_tot, n_suc, conf, digits = 2) {
  
  stopifnot(n_tot > 0 && n_suc > 0 && n_tot >= n_suc)
  stopifnot(conf > 0 && conf < 1)
  
  bsc_ci <- .bscCI(n_tot, n_suc, conf)
  if (!missing(digits)) {
    round(bsc_ci, digits)
  }
  else {
    return(bsc_ci)
  }
}
