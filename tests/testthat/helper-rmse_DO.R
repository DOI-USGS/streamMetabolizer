#' Calculate the RMSE of DO predictions relative to observations
#' 
#' @param DO_preds predictions from predict_DO()
rmse_DO <- function(DO_preds) {
  sqrt(mean((DO_preds$DO.obs - DO_preds$DO.mod)^2, na.rm=TRUE))
}