#' Determine how many cores to use for an MCMC run
#'
#' Shared by \code{runstan_bayes} (one-station) and \code{metab_bayes_2s}
#' (two-station): detects the number of cores available on the machine,
#' falls back to 1 if detection fails, and caps the requested core count at
#' whatever is actually available.
#'
#' @param n_cores the number of cores requested for this run
#' @param n_chains the number of chains being requested, used only to
#'   reconstruct the verbose status message; if NULL (the default), the
#'   chains clause is omitted from the message. Ignored when
#'   \code{verbose=FALSE}
#' @param verbose logical. if TRUE, emit a status message reporting the
#'   number of cores requested vs. available
#' @return the number of cores to actually use, i.e.
#'   \code{min(detected_cores, n_cores)}
#' @keywords internal
mm_determine_cores <- function(n_cores, n_chains=NULL, verbose=FALSE) {
  tot_cores <- parallel::detectCores()
  if(!is.finite(tot_cores)) tot_cores <- 1
  n_cores <- min(tot_cores, n_cores)
  if(verbose) {
    chains_clause <- if(is.null(n_chains)) "" else paste0(n_chains," chains on ")
    message(paste0("MCMC (","Stan","): requesting ",chains_clause,n_cores," of ",tot_cores," available cores"))
  }
  n_cores
}
