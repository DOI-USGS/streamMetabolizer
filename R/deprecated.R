#' Warn about deprecation of a units-related argument or function
#' @keywords internal
unitted_deprecate_warn <- function(what) {
  deprecate_warn(
    what = sprintf('streamMetabolizer::%s', what),
    when = "0.12.0",
    details = "streamMetabolizer will someday stop using units and the unitted package.")
}
