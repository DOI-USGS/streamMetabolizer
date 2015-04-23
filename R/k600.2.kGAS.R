#' Returns the gas exchange velocity for gas of interest w/ no unit conversions
#' 
#' @importFrom LakeMetabolizer k600.2.kGAS
k600.2.kGAS = function(ts.data, gas="O2"){
  LakeMetabolizer::k600.2.kGAS(ts.data, gas)
}