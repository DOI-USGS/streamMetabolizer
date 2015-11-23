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
plot_DO_preds <- function(DO_preds, plot_as=c('conc','pctsat','ddodt')) {
  
  plot_as <- match.arg(plot_as, several.ok=TRUE)
  
  DO.obs <- DO.mod <- DO.sat <- '.dplyr.var'
  DO_preds_conc <- mutate(
    DO_preds, as='conc', var='DO (mg/L)', col='red4', lab='DO (mg/L)',
    mod=DO.mod, 
    obs=DO.obs)
  DO_preds_pctsat <- mutate(
    DO_preds, as='pctsat', var='DO (% sat)', col='forestgreen', lab='DO (% sat)',
    mod=100*DO.mod/DO.sat, 
    obs=100*DO.obs/DO.sat)
  DO_preds_ddodt <- 
    mutate(
      DO_preds[-1,], as='ddodt', var='dDO/dt (mg/L/d)', col='navy', lab='dDO/dt~(mg~L^-1~d^-1)',
      mod = diff(DO_preds$DO.mod)/as.numeric(diff(DO_preds$local.time), units="days"),
      obs = diff(DO_preds$DO.obs)/as.numeric(diff(DO_preds$local.time), units="days")) %>%
    mutate(
      mod = ifelse(diff(DO_preds$local.date)==0, mod, NA),
      obs = ifelse(diff(DO_preds$local.date)==0, obs, NA))
  
  DO_preds_all <- bind_rows(DO_preds_conc, DO_preds_pctsat, DO_preds_ddodt) %>%
    mutate(var=ordered(var, c(conc='DO (mg/L)', pctsat='DO (% sat)', ddodt='dDO/dt (mg/L/d)')[plot_as]))
  
  local.time <- mod <- local.date <- col <- obs <- '.ggplot.var'
  g <- ggplot(v(DO_preds_all)[DO_preds_all$as %in% plot_as,], aes(x=local.time, group=local.date, color=col)) +
    geom_line(aes(y=mod), size=0.8) +
    geom_point(aes(y=obs), alpha=0.3) +
    scale_color_identity(guide='none') +
    theme_bw() + facet_grid(var ~ ., scales="free_y") + 
    xlab('Local time') + ylab("Predictions (lines) and observations (points)")
  
  suppressWarnings(print(g))
}