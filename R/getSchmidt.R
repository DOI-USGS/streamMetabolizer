#' Returns Schmidt number for a specific gas at a given temperature
#' 
#' @importFrom LakeMetabolizer getSchmidt
getSchmidt  <-	function(temperature, gas){
  LakeMetabolizer::getSchmidt(temperature, gas)
}