#' Build a closure that pivots a per-row vector into a time-by-date matrix
#'
#' Both \code{prepdata_bayes} (one-station) and \code{prepdata_bayes_2s}
#' (two-station) reshape several per-row vectors (DO, depth, light, etc.)
#' into an obs-per-day x num-days matrix, assuming that the input vector is
#' already sorted so that each date occupies a contiguous block of
#' \code{n_per_group} rows. This function returns the reshaping closure;
#' \code{\link{mm_check_dates_contiguous}} verifies the contiguous-block
#' assumption actually holds.
#'
#' @param n_per_group the number of rows per date (must be the same for
#'   every date; callers are responsible for having already confirmed this)
#' @param n_groups the number of distinct dates
#' @return a function of one argument, \code{vec}, that reshapes \code{vec}
#'   into an \code{n_per_group} x \code{n_groups} matrix, filling by column
#' @keywords internal
mm_time_by_date_matrix <- function(n_per_group, n_groups) {
  function(vec) matrix(vec, nrow=n_per_group, ncol=n_groups, byrow=FALSE)
}

#' Confirm that each date occupies a contiguous block of rows
#'
#' Used immediately after pivoting a date vector with the closure from
#' \code{\link{mm_time_by_date_matrix}}, to catch input that wasn't actually
#' sorted by date/time before pivoting (which would otherwise let the matrix
#' reshape silently scramble which rows belong to which date). Shared by
#' \code{prepdata_bayes} and \code{prepdata_bayes_2s}.
#'
#' @param date_mat the date-identifier vector (as used to build
#'   \code{date_table}) already pivoted via
#'   \code{\link{mm_time_by_date_matrix}}'s closure
#' @param date_table a table of date counts, as from \code{table(date_vec)},
#'   giving the expected date for each column of \code{date_mat}
#' @param error_message the message to pass to \code{stop()} if the dates
#'   are not contiguous; left to the caller so each can keep its own wording
#' @keywords internal
mm_check_dates_contiguous <- function(date_mat, date_table, error_message) {
  unique_dates_per_col <- apply(date_mat, MARGIN=2, FUN=unique)
  if(is.list(unique_dates_per_col) || !isTRUE(all.equal(unname(unique_dates_per_col), names(date_table)))) {
    stop(error_message)
  }
  invisible(TRUE)
}
