#' Nighttime regression estimation of ER and K600
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_night_basic <- function(
  # possible additional args: 
  # filter_width=3 (in night_dat$DO.obs.smooth <- u(c(stats::filter(night_dat$DO.obs, rep(1/3, 3), sides=2)), get_units(night_dat$DO.obs)))
  # night_end or night_length (time or time period after which to cut off nighttime data. NA could mean use full night to dawn)
) {
  
  as.list(environment())
  
}