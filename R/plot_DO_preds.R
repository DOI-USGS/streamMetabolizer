#' Plot predictions produced with predict_DO
#' 
#' Plots modeled values as lines, observed values as points
#' 
#' @param DO_preds a data.frame of predictions such as that returned by 
#'   predict_DO()
#' @param plot_as character. Should the plot display predicted & observed values
#'   in concentration (conc) or as percent of saturation (pctsat)? The default
#'   is to plot both.
#' @examples 
#' \dontrun{
#' mm <- metab_night(v(french))
#' plot_DO_preds(predict_DO(mm)[1:360,])
#' }
#' @import ggplot2
#' @import dplyr
#' @export
plot_DO_preds <- function(DO_preds, plot_as=c('both','conc','pctsat')) {
  
  plot_as <- match.arg(plot_as)
  
  DO.mod <- '.ggplot.var'
  g <- switch(
    plot_as,
    'conc'={
      ggplot(v(DO_preds), aes(x=local.time)) + 
        geom_line(aes(y=DO.mod, group=local.date), color='maroon', size=0.8) +
        geom_point(aes(y=DO.obs), color='navy', alpha=0.5) +
        theme_bw() + xlab('Local time') + ylab('DO (mg/L)') 
    },
    'pctsat'={
      ggplot(v(DO_preds), aes(x=local.time)) + 
        geom_line(aes(y=100*DO.mod/DO.sat, group=local.date), color='purple3', size=0.8) +
        geom_point(aes(y=100*DO.obs/DO.sat), color='forestgreen', alpha=0.5) +
        theme_bw() + xlab('Local time') + ylab('DO (% sat)')
    },
    'both'={
      DO_recast <- rbind(
        mutate(DO_preds, var='DO (mg/L)', DO.mod=DO.mod, DO.obs=DO.obs),
        mutate(DO_preds, var='DO (% sat)', DO.mod=100*DO.mod/DO.sat, DO.obs=100*DO.obs/DO.sat))
      ggplot(DO_recast, aes(x=local.time, y=DO.mod, group=local.date)) + 
        geom_line(aes(color=var), size=0.8) +
        geom_point(aes(y=DO.obs, color=var), alpha=0.5) +
        scale_color_discrete(guide='none') +
        theme_bw() + facet_grid(var ~ ., scales="free_y") + 
        xlab('Local time') + ylab("DO predictions (% sat above, mg/L below)")
    })
  suppressWarnings(print(g))
}