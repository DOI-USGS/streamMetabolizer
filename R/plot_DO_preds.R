#' Plot predictions produced with predict_DO
#' 
#' Plots modeled values as lines, observed values as points
#' 
#' @param DO_preds a data.frame of predictions such as that returned by 
#'   predict_DO()
#' @param y_var character. Should the plot display predicted & observed values 
#'   in concentration (conc) or as percent of saturation (pctsat)? The default 
#'   is to plot both.
#' @param style character indicating which graphics package to use
#' @param y_lim list of named vectors, each of which has length 2 and is numeric
#'   and has a name in the possible values of y_var. NA within a vector
#'   indicates that the data range should be used. for ggplot2, y_lim is only
#'   used to exclude values outside that range and is ignored if the data span a
#'   narrower range
#' @examples 
#' \dontrun{
#' mm <- metab_night(specs(mm_name('night')), data=data_metab('3', day_start=12, day_end=36))
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(predict_DO(mm), style='dygraphs', y_var='pctsat')
#' }
#' @import dplyr
#' @importFrom unitted v
#' @export
plot_DO_preds <- function(DO_preds, y_var=c('conc','pctsat','ddodt'), 
                          style=c('ggplot2','dygraphs'),
                          y_lim=list(conc=c(NA,NA), pctsat=c(NA,NA), ddodt=c(NA,NA))) {
  
  style <- match.arg(style)
  y_var <- match.arg(y_var, several.ok=TRUE)
  
  params <- list(
    xlab='Local time',
    ylab='Predictions (lines) and observations (points)',
    colors=list(conc=c('#A64B00','#FF7400'), pctsat=c('#007929','#23BC47'), ddodt=c('#05326D','#4282D3'))
  )
  
  DO.obs <- DO.mod <- DO.sat <- '.dplyr.var'
  DO_preds_conc <- mutate(
    DO_preds, as='conc', var='DO (mg/L)', col1=params$colors$conc[1], col2=params$colors$conc[2], lab='DO (mg/L)',
    mod=DO.mod, 
    obs=DO.obs)
  DO_preds_pctsat <- mutate(
    DO_preds, as='pctsat', var='DO (% sat)', col1=params$colors$pctsat[1], col2=params$colors$pctsat[2], lab='DO (% sat)',
    mod=100*DO.mod/DO.sat, 
    obs=100*DO.obs/DO.sat)
  DO_preds_ddodt <- 
    mutate(
      DO_preds[-1,], as='ddodt', var='dDO/dt (mg/L/d)', col1=params$colors$ddodt[1], col2=params$colors$ddodt[2], lab='dDO/dt~(mg~L^-1~d^-1)',
      mod = diff(DO_preds$DO.mod)/as.numeric(diff(DO_preds$solar.time), units="days"),
      obs = diff(DO_preds$DO.obs)/as.numeric(diff(DO_preds$solar.time), units="days")) %>%
    mutate(
      mod = ifelse(diff(DO_preds$date)==0, mod, NA),
      obs = ifelse(diff(DO_preds$date)==0, obs, NA))
  
  DO_preds_all <- bind_rows(DO_preds_conc, DO_preds_pctsat, DO_preds_ddodt) %>%
    mutate(var=ordered(var, c(conc='DO (mg/L)', pctsat='DO (% sat)', ddodt='dDO/dt (mg/L/d)')[y_var]))
  
  plot_out <- switch(
    style,
    'ggplot2' = {
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      
      . <- solar.time <- mod <- date <- col1 <- col2 <- obs <- '.ggplot.var'
      preds_ggplot <- v(DO_preds_all) %>%
        filter(as %in% y_var)
      if('conc' %in% names(y_lim)) {
        lim <- y_lim[['conc']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'conc' | (mod >= lim & obs >= lim))
        lim <- y_lim[['conc']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'conc' | (mod <= lim & obs <= lim))
      }
      if('pctsat' %in% names(y_lim)) {
        lim <- y_lim[['pctsat']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'pctsat' | (mod >= lim & obs >= lim))
        lim <- y_lim[['pctsat']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'pctsat' | (mod <= lim & obs <= lim))
      }
      if('ddodt' %in% names(y_lim)) {
        lim <- y_lim[['ddodt']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ddodt' | (mod >= lim & obs >= lim))
        lim <- y_lim[['ddodt']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ddodt' | (mod <= lim & obs <= lim))
      }
      g <- ggplot2::ggplot(preds_ggplot, ggplot2::aes(x=solar.time, group=date)) +
        ggplot2::geom_point(ggplot2::aes(y=obs, color=col2), alpha=0.6) +
        ggplot2::geom_line(ggplot2::aes(y=mod, color=col1), size=0.8) +
        ggplot2::scale_color_identity(guide='none') +
        ggplot2::theme_bw() + 
        ggplot2::facet_grid(var ~ ., scales="free_y") + 
        ggplot2::xlab(params$xlab) + ggplot2::ylab(params$ylab)
      
      suppressWarnings(print(g))
    },
    'dygraphs' = {
      if(!requireNamespace("dygraphs", quietly=TRUE))
        stop("call install.packages(dygraphs') before plotting with style='dygraphs'")
      if(!requireNamespace("xts", quietly=TRUE))
        stop("call install.packages(dygraphs') before plotting with style='dygraphs'")
      
      . <- '.dplyr.var'
      preds_xts <- v(DO_preds_all) %>%
        filter(as %in% y_var) %>%
        arrange(solar.time) %>%
        group_by(date) %>%
        do(., {
          out <- .[c(seq_len(nrow(.)),nrow(.)),]
          out[nrow(.)+1,c('mod','obs')] <- NA
          out
        }) %>%
        ungroup()
      
      prep_dygraph <- function(y_var) { 
        preds_xts %>% 
          filter(as==y_var) %>% 
          select(mod,obs,solar.time) %>%
          mutate(solar.time=lubridate::force_tz(solar.time, Sys.getenv("TZ"))) %>% # dygraphs makes some funky tz assumptions. this seems to help.
          xts::xts(x=select(., -solar.time), order.by=.$solar.time, unique=FALSE, tzone=Sys.getenv("TZ"))
      }
      if(length(y_var) > 1) warning("can only plot one dygraph y_var at a time for now; plotting y_vars in succession")
      sapply(y_var, function(yvar) {
        y_var_long <- preds_xts %>% filter(as==yvar) %>% slice(1) %>% .[['var']] %>% as.character()
        y_var_col <- params$colors[[yvar]]
        dat <- prep_dygraph(yvar)
        ymin <- max(c(min(c(dat[,1], dat[,2]), na.rm=TRUE), y_lim[[yvar]][1]), na.rm=TRUE)
        ymax <- min(c(max(c(dat[,1], dat[,2]), na.rm=TRUE), y_lim[[yvar]][2]), na.rm=TRUE)
        dygraphs::dygraph(dat, xlab=params$xlab, ylab=y_var_long, group='plot_DO_preds') %>%
          dygraphs::dySeries('mod', drawPoints = FALSE, label=paste0("Modeled ", y_var_long), color=y_var_col[1]) %>%
          dygraphs::dySeries('obs', drawPoints = TRUE, strokeWidth=0, label=paste0("Observed ", y_var_long), color=y_var_col[2]) %>%
          dygraphs::dyAxis('y', valueRange=(c(ymin,ymax)+(ymax-ymin)*c(-0.05,0.15))) %>%
          dygraphs::dyOptions(colorSaturation=1) %>%
          dygraphs::dyLegend(labelsSeparateLines = TRUE, width=300) %>%
          dygraphs::dyRangeSelector(height = 20) %>%
          print()
      })
    }
  )
  
  invisible(plot_out)
}