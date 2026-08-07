#' Placeholder daily estimates for dates that produced no fit
#'
#' The nine GPP/ER/K600 credible-interval columns that
#' \code{\link{predict_metab}} and \code{\link{get_params}} look for, filled
#' with \code{NA} for each of \code{dates}, so that a failed fit still
#' reports which dates were attempted instead of returning nothing.
#'
#' Shared by \code{bayes_allply} (one-station, all dates in one Stan call)
#' and the two-station fitting functions (see \code{metab_bayes_2s.R}), which
#' built byte-identical frames independently before this was extracted.
#' Keeping the column list in one place is the point: it also appears in
#' \code{bayes_1ply} and in \code{\link{metab_bayes}}'s removed-dates join,
#' and those two are deliberately \emph{not} routed through here -- the
#' former returns a single row with no \code{date} column (its caller,
#' \code{\link{mm_model_by_ply}}, attaches one), and the latter builds its
#' columns onto an existing data.frame of removed dates and marks them
#' \code{valid_day=FALSE}. Serving either from here would take more
#' branching than the duplication costs.
#'
#' @param dates Date vector of the dates to report as attempted
#' @return data.frame with one row per date: \code{date} plus the 2.5/50/97.5
#'   percentile columns for \code{GPP_daily}, \code{ER_daily}, and
#'   \code{K600_daily}, all \code{NA}
#' @keywords internal
mm_na_daily <- function(dates) {
  na_vec <- rep(as.numeric(NA), length(dates))
  data.frame(
    date=dates,
    GPP_daily_2.5pct=na_vec, GPP_daily_50pct=na_vec, GPP_daily_97.5pct=na_vec,
    ER_daily_2.5pct=na_vec, ER_daily_50pct=na_vec, ER_daily_97.5pct=na_vec,
    K600_daily_2.5pct=na_vec, K600_daily_50pct=na_vec, K600_daily_97.5pct=na_vec)
}
