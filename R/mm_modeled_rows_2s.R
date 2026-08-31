#' @include mm_lag_2s.R
NULL

#' Resolve a two-station alignment into the rows that will be modeled
#'
#' Produces one row per modeled observation, each column drawn from the rows
#' that actually supply it: upstream DO and DO saturation from
#' \code{shift_idx} (the rows one travel time earlier), everything else from
#' \code{keep}. This is the one definition of the data a two-station fit sees,
#' so that everything acting on those values -- reshaping them for the model,
#' testing them for validity -- agrees on which value belongs to which row.
#'
#' Columns keep their input names rather than the model's internal ones, so
#' that anything reporting on this frame names a column the user supplied
#' (\code{DO.obs.up}, not \code{DO_obs_up}).
#'
#' @param data data.frame as validated by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_bayes_2s}}, units optional. Always the \emph{full}
#'   dataset, never a per-day slice: \code{aln}'s indices point into it, and
#'   \code{shift_idx} routinely reaches into earlier days' rows.
#' @param aln an alignment as returned by \code{mm_align_2s}, or a single-day
#'   slice of one. Supplying an alignment that doesn't correspond to
#'   \code{data} will silently produce wrong rows.
#' @return a data.frame, units stripped, carrying the two-station data columns
#'   with one row per modeled observation, parallel to \code{aln$keep} and
#'   \code{aln$date}
#' @importFrom unitted v
#' @keywords internal
mm_modeled_rows_2s <- function(data, aln) {

  # units are stripped here rather than by each caller: Stan cannot take
  # unitted matrices, and the validity tests are simpler on plain numerics
  data <- v(data)

  keep <- aln$keep
  shift_idx <- aln$shift_idx

  data.frame(
    solar.time  = data$solar.time[keep],
    DO.obs.up   = data$DO.obs.up[shift_idx],
    DO.sat.up   = data$DO.sat.up[shift_idx],
    DO.obs.down = data$DO.obs.down[keep],
    DO.sat.down = data$DO.sat.down[keep],
    light       = data$light[keep],
    depth       = data$depth[keep],
    temp.water  = data$temp.water[keep],
    travel.time = data$travel.time[keep]
  )
}
