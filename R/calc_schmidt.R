#' @title Returns Schmidt number for a specific gas at a given temperature
#' @param temp.water water temperature in degC
#' @param gas string for gas code. Valid inputs include: He, O2, CO2, CH4, SF6, N2O, Ar, and N2
#' @return Schmidt number (unitless)
#' @note
#' Temperature range is only valid from 4-35 deg Celsius
#' @references
#' Raymond, Peter A., Christopher J. Zappa, David Butman, Thomas L. Bott, Jody Potter, 
#' Patrick Mulholland, Andrew E. Laursen, William H. McDowell, and Denis Newbold. 
#' Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers. 
#' Limnology & Oceanography: Fluids & Environments 2 (2012): 41-53.
#' @examples
#' calc_schmidt(temp.water = 25, gas = 'O2')
#' temp.water <- unitted::u(25, 'degC')
#' calc_schmidt(temp.water, gas = 'O2')
#' @importFrom LakeMetabolizer getSchmidt
#' @importFrom unitted u verify_units is.unitted
#' @export
calc_schmidt  <-	function(temp.water, gas){
  
  schmidt <- LakeMetabolizer::getSchmidt(temperature = v(temp.water), gas)
  
  if (is.unitted(temp.water) && verify_units(temp.water,'degC', list(TRUE,FALSE))){
    return(u(schmidt, ""))
  } else {
    return(schmidt)
  }

}