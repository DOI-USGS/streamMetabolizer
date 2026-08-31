#' Attach fitting warnings and errors to a predictions data.frame
#'
#' Adds the date-specific fitting messages as \code{msgs.fit}, then the
#' prediction-level \code{warnings} and \code{errors} columns. General errors
#' during fitting almost always prohibit prediction and general warnings
#' almost always don't, so both columns are keyed off the overall errors:
#' predictions are NA where fitting errored and unremarkable otherwise. There
#' is no attempt to separate prediction-specific messages from fitting ones,
#' since both come out of the same engine and arrive indistinguishable.
#'
#' Shared by the Bayesian model types' \code{\link{predict_metab}} methods,
#' which do not call \code{\link{get_params}} and so must attach the overall
#' messages themselves.
#'
#' @param preds data.frame of predictions, with a \code{date} column.
#' @param fit the fit's per-date data.frame, already filtered to the dates
#'   being predicted.
#' @param warnings.overall character vector of the fit's overall warnings.
#' @param errors.overall character vector of the fit's overall errors.
#' @return \code{preds} with \code{msgs.fit}, \code{warnings}, and
#'   \code{errors} columns added.
#' @keywords internal
mm_attach_fit_msgs <- function(preds, fit, warnings.overall, errors.overall) {

  warnings <- errors <- '.dplyr.var'
  if(!is.null(fit) && all(c('date','warnings','errors') %in% names(fit))) {
    messages <- fit %>%
      select(date, warnings, errors) %>%
      compress_msgs('msgs.fit', warnings.overall=warnings.overall, errors.overall=errors.overall)
    preds <- full_join(preds, messages, by='date', copy=TRUE)
  } else {
    preds <- mutate(preds, msgs.fit=NA)
  }

  mutate(
    preds,
    warnings=if(length(errors.overall) > 0) NA else '',
    errors=if(length(errors.overall) > 0) NA else '')
}
