#' Check the format of keep_mcmcs and keep_mcmc_data
#'
#' Either option may be a single logical, or a vector of dates naming the
#' fits to keep. The date-vector form only means something when each date has
#' its own fit, so it is rejected unless \code{split_dates} is TRUE; where it
#' is allowed, it is coerced to \code{Date}.
#'
#' Shared by the Bayesian model types so that they cannot disagree about what
#' a legal value is.
#'
#' @param specs a specs list, as from \code{\link{specs}}.
#' @return \code{specs}, with \code{keep_mcmcs} and \code{keep_mcmc_data}
#'   coerced to \code{Date} if they were given as date vectors.
#' @keywords internal
mm_check_keep_mcmc_specs <- function(specs) {

  if(is.logical(specs$keep_mcmcs)) {
    if(length(specs$keep_mcmcs) != 1) {
      stop("if keep_mcmcs is logical, it must have length 1")
    }
  } else if(specs$split_dates == FALSE) {
    stop("if split_dates==FALSE, keep_mcmcs must be a single logical value")
  }
  if(is.logical(specs$keep_mcmc_data)) {
    if(length(specs$keep_mcmc_data) != 1) {
      stop("if keep_mcmc_data is logical, it must have length 1")
    }
  } else if(specs$split_dates == FALSE) {
    stop("if split_dates==FALSE, keep_mcmc_data must be a single logical value")
  }

  # only reachable when split_dates is TRUE: a non-logical value has already
  # been rejected above otherwise
  if(!is.logical(specs$keep_mcmcs)) specs$keep_mcmcs <- as.Date(specs$keep_mcmcs)
  if(!is.logical(specs$keep_mcmc_data)) specs$keep_mcmc_data <- as.Date(specs$keep_mcmc_data)

  specs
}
