#' Calculates the equilibrium saturation concentration of oxygen in water at the supplied conditions
#' 
#' @importFrom LakeMetabolizer o2.at.sat
o2.at.sat <- function(ts.data, baro, altitude=0, salinity=0, model='garcia'){
  LakeMetabolizer::o2.at.sat(ts.data, baro, altitude, salinity, model)
}