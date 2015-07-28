#' Plot predictions produced with predict_DO
#' 
#' Plots modeled values as lines, observed values as points
#' 
#' @param DO_preds a data.frame of predictions such as that returned by
#'   predict_DO()
#' @param plot_as character. Should the plot display predicted & observed values
#'   in concentration (conc) or as percent of saturation (pctsat)?
#' @examples 
#' mm <- metab_night_data=v(french))
#' plot_DO_preds(predict_DO(mm)[1:360,])
#' @export
plot_DO_preds <- function(DO_preds, plot_as=c('conc','pctsat')) {
  
  plot_as <- match.arg(plot_as)
  
  if(plot_as=='conc') {
    g <- ggplot(v(DO_preds), aes(x=local.time)) + 
      geom_line(aes(y=DO.mod, group=date), color='maroon', size=0.8) +
      geom_point(aes(y=DO.obs), color='navy', alpha=0.5) +
      theme_bw() + ylab('DO (mg/L)')
  
  } else if(plot_as=='pctsat') {
    g <- ggplot(v(DO_preds), aes(x=local.time)) + 
      geom_line(aes(y=100*DO.mod/DO.sat, group=date), color='purple3', size=0.8) +
      geom_point(aes(y=100*DO.obs/DO.sat), color='forestgreen', alpha=0.5) +
      theme_bw() + ylab('DO (% sat)')
  }
  print(g)
}