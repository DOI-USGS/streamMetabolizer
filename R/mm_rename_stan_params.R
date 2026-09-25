#' Rename Stan output columns to streamMetabolizer conventions
#'
#' Stan prohibits '.' in variable names, so parameters come back with '_'
#' where the rest of the package uses '.'. Converts the parameter names
#' themselves, then the summary-statistic suffixes.
#'
#' Shared by the Bayesian model types' \code{\link{get_params}} methods.
#'
#' @param df data.frame of per-date Stan output, as in a fit's \code{daily}.
#' @param params_out character vector of the Stan parameter names to convert,
#'   as in \code{specs$params_out}.
#' @return \code{df}, with its columns renamed.
#' @keywords internal
mm_rename_stan_params <- function(df, params_out) {

  parnames <- setNames(gsub('_', '\\.', params_out), params_out)
  # longest first, so a shorter name can never rewrite part of a longer one
  parnames <- parnames[order(nchar(parnames), decreasing=TRUE)]
  for(i in seq_along(parnames)) {
    names(df) <- gsub(names(parnames[i]), parnames[[i]], names(df))
  }
  names(df) <- gsub('_mean$', '', names(df))
  names(df) <- gsub('_sd$', '.sd', names(df))
  names(df) <- gsub('_50pct$', '.median', names(df))
  names(df) <- gsub('_2.5pct$', '.lower', names(df))
  names(df) <- gsub('_97.5pct$', '.upper', names(df))

  df
}
