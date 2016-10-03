#' Deprecated Functions in package streamMetabolizer
#' 
#' These functions are provided for compatibility with older versions of 
#' \code{streamMetabolizer} only, and may be defunct as soon as the next 
#' release.
#' 
#' \itemize{
#'   \item \code{\link{calc_DO_deficit}} - instead, subtract \code{DO.obs} from output of \code{\link{calc_DO_sat}}
#'   \item \code{calc_DO_at_sat} - use \code{\link{calc_DO_sat}} instead
#'   \item \code{\link{calc_is_daytime}} - if you like and want this function, submit a GitHub issue to keep it
#'   \item \code{\link{calc_sun_rise_set}} - if you like and want this function, submit a GitHub issue to keep it
#'   \item \code{\link{calc_velocity}} - if you like and want this function, submit a GitHub issue to keep it
#' }
#' 
#' These functions will are currently exported but will soon be internal-only.
#' You are encouraged not to use these.
#' 
#' \itemize{
#'   \item \code{\link{convert_date_to_doyhr}}
#'   \item \code{\link{convert_doyhr_to_date}}
#'   \item \code{\link{lookup_google_timezone}} - use \code{\link{lookup_timezone}} instead
#' }
#' 
#' @name streamMetabolizer-deprecated
NULL
