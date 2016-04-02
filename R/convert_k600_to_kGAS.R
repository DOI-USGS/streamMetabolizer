#' Returns the gas exchange velocity for gas of interest w/ no unit conversions
#' 
#' @param k600 k600 as vector of numbers or single number
#' @param temperature Water temperature (deg C) as vector array of numbers or single number
#' @param gas gas for conversion, as string (e.g., 'CO2' or 'O2')
#' @return Numeric value of gas exchange velocity for gas
#' 
#' @importFrom LakeMetabolizer k600.2.kGAS.base
#' @importFrom unitted is.unitted verify_units v u get_units
#' @export
convert_k600_to_kGAS = function(k600, temperature, gas="O2") {
  # units checks
  with_units <- is.unitted(k600) || is.unitted(temperature)
  if(with_units) {
    if(!is.unitted(k600)) stop("if temperature is unitted, then k600 should also be unitted")
    if(!is.unitted(temperature)) stop("if k600 is unitted, then temperature should also be unitted")
    verify_units(temperature, 'degC') 
  }
  
  # suppressing "In getSchmidt(temperature, gas) : temperature out of range" b/c it's way too common
  out <- suppressWarnings(LakeMetabolizer::k600.2.kGAS.base(v(k600), v(temperature), v(gas)))
  
  if(with_units) u(out, get_units(k600)) else out
}

#' Returns the gas exchange velocity as k600 for gas of interest w/ no unit conversions
#' 
#' @param kGAS k of gas as vector of numbers or single number
#' @param temperature Water temperature (deg C) as vector array of numbers or single number
#' @param gas gas for conversion, as string (e.g., 'CO2' or 'O2')
#' @return Numeric value of gas exchange velocity for gas
#' 
#' @importFrom LakeMetabolizer k600.2.kGAS.base
#' @importFrom unitted is.unitted verify_units v
#' @export
convert_kGAS_to_k600 = function(kGAS, temperature, gas="O2") {
  # units checks
  with_units <- is.unitted(kGAS) || is.unitted(temperature)
  if(with_units) {
    if(!is.unitted(kGAS)) stop("if temperature is unitted, then kGAS should also be unitted")
    if(!is.unitted(temperature)) stop("if kGAS is unitted, then temperature should also be unitted")
    verify_units(temperature, 'degC') 
  }
  
  # suppressing "In getSchmidt(temperature, gas) : temperature out of range" b/c it's way too common
  conversion <- 1/suppressWarnings(LakeMetabolizer::k600.2.kGAS.base(1, v(temperature), v(gas)))
  kGAS * conversion # do the conversion and adopt the units of kGAS, if any
}
