#' Convert from photosynthetically active to shortwave radiation
#' 
#' Convert photosynthetically active radiation (PAR) to shortwave radiation 
#' (SW).
#' 
#' @param par Vector of photosynthetically active radiation (400-700 nm;
#'   umol/m^2/sec)
#' @param coeff Numerical coefficient to convert PAR (umol/m^2/sec) to SW
#'   (W/m^2). Defaults to value from Britton and Dodd (1976).
#' @return Numeric vector of shortwave values with units W/m^2
#'
#' @examples
#' convert_PAR_to_SW(par=400, coef=0.47)
#' 
#' @importFrom LakeMetabolizer par.to.sw.base
#' @export
convert_PAR_to_SW <- function(par, coeff=0.473) {
  (LakeMetabolizer::par.to.sw.base(par, coeff)) # parentheses to return value as visible
}

#' Convert from shortwave to photosynthetically active radiation
#' 
#' Convert shortwave radiation (SW) to photosynthetically active radiation
#' (PAR).
#' 
#' @param sw Vector of shortwave radiation (W/m^2)
#' @param coeff Numerical coefficient to convert SW (W/m^2) to PAR 
#'   (umol/m^2/sec). Defaults to value from Britton and Dodd (1976).
#' @return Numeric vector of PAR values in units umol/m^2/sec
#' 
#' @examples
#' convert_SW_to_PAR(sw=800)
#' convert_SW_to_PAR(sw=800, coef=2.1)
#' 
#' @importFrom LakeMetabolizer sw.to.par.base
#' @export
convert_SW_to_PAR <- function(sw, coeff=2.114) {
  (LakeMetabolizer::sw.to.par.base(sw, coeff)) # parentheses to return value as visible
}
