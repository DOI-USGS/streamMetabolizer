#' \code{specs_night_basic} - nighttime regression estimation of ER and K600
#' 
#' @rdname specs_night
#'   
#' @param calc_DO_fun The function to use in calculating DO; probably either 
#'   \code{calc_DO_mod} or \code{calc_DO_mod_by_diff}. For nighttime regression,
#'   this function is not used for model fitting at all, but is still used to 
#'   re-predict the trajectory of DO concentrations from the fitted values of 
#'   K600 and ER. It is unclear to the package authors whether calc_DO_mod or 
#'   calc_DO_mod_by_diff is more appropriate.
#' @export
#' @family model_specs
specs_night_basic <- function(
  calc_DO_fun=calc_DO_mod
  # possible additional args: 
  # filter_width=3 (in night_dat$DO.obs.smooth <- u(c(stats::filter(night_dat$DO.obs, rep(1/3, 3), sides=2)), get_units(night_dat$DO.obs)))
  # night_end or night_length (time or time period after which to cut off nighttime data. NA could mean use full night to dawn)
) {
  list(
    calc_DO_fun=calc_DO_fun
  )
}