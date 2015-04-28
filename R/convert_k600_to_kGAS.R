#' Returns the gas exchange velocity for gas of interest w/ no unit conversions
#' 
#' @param k600 k600 as vector of numbers or single number
#' @param temperature Water temperature (deg C) as vector array of numbers or single number
#' @param gas gas for conversion, as string (e.g., 'CO2' or 'O2')
#' @return Numeric value of gas exchange velocity for gas
#' 
#' @importFrom LakeMetabolizer k600.2.kGAS.base
#' @export
convert_k600_to_kGAS = function(k600, temperature, gas="O2") {
  # units checks would be good here
  LakeMetabolizer::k600.2.kGAS.base(k600, temperature, gas)
}

#' Returns the gas exchange velocity as k600 for gas of interest w/ no unit conversions
#' 
#' @param kGAS k of gas as vector of numbers or single number
#' @param temperature Water temperature (deg C) as vector array of numbers or single number
#' @param gas gas for conversion, as string (e.g., 'CO2' or 'O2')
#' @return Numeric value of gas exchange velocity for gas
#' 
#' @importFrom LakeMetabolizer k600.2.kGAS.base
#' @export
convert_kGAS_to_k600 = function(kGAS, temperature, gas="O2") {
  # units checks would be good here
  conversion <- 1/LakeMetabolizer::k600.2.kGAS.base(1, temperature, gas)
  kGAS * conversion
}
