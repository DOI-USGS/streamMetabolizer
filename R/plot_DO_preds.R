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
  
  if(is(DO_preds, 'metab_model')) DO_preds <- predict_DO(DO_preds)
  
  style <- match.arg(style)
  y_var <- match.arg(y_var, several.ok=TRUE)
  
  params <- list(
    xlab='Local time',
    ylab='Predictions (lines) and observations (points)',
    colors=list(conc=c('#CE9C59', '#A64B00','#FF7400'), 
                pctsat=c('#7CA586','#007929','#23BC47'), 
                ddodt=c('#4A5869','#05326D','#4282D3'))
  )
  
  DO.obs <- DO.mod <- DO.sat <- '.dplyr.var'
  DO_preds_conc <- mutate(
    DO_preds, as='conc', var='DO (mg/L)', lab='DO (mg/L)',
    col.pure=params$colors$conc[1], col.mod=params$colors$conc[2], col.obs=params$colors$conc[3], 
    pure=if(exists('DO.pure', DO_preds)) DO.pure else NA,
    mod=DO.mod, 
    obs=DO.obs)
  DO_preds_pctsat <- mutate(
    DO_preds, as='pctsat', var='DO (% sat)', lab='DO (% sat)',
    col.pure=params$colors$pctsat[1], col.mod=params$colors$pctsat[2], col.obs=params$colors$pctsat[3], 
    pure=if(exists('DO.pure', DO_preds)) 100*DO.pure/DO.sat else NA, 
    mod=100*DO.mod/DO.sat, 
    obs=100*DO.obs/DO.sat)
  DO_preds_ddodt <- 
    mutate(
      DO_preds[-1,], as='ddodt', var='dDO/dt (mg/L/d)', lab='dDO/dt~(mg~L^-1~d^-1)',
      col.pure=params$colors$ddodt[1], col.mod=params$colors$ddodt[2], col.obs=params$colors$ddodt[3], 
      pure = if(exists('DO.pure', DO_preds)) diff(DO_preds$DO.pure)/as.numeric(diff(DO_preds$solar.time), units="days") else NA,
      mod = diff(DO_preds$DO.mod)/as.numeric(diff(DO_preds$solar.time), units="days"),
      obs = diff(DO_preds$DO.obs)/as.numeric(diff(DO_preds$solar.time), units="days")) %>%
    mutate(
      pure = ifelse(diff(DO_preds$date)==0, pure, NA),
      mod = ifelse(diff(DO_preds$date)==0, mod, NA),
      obs = ifelse(diff(DO_preds$date)==0, obs, NA))
  
  var <- '.dplyr.var'
  DO_preds_all <- bind_rows(DO_preds_conc, DO_preds_pctsat, DO_preds_ddodt) %>%
    mutate(var=ordered(var, c(conc='DO (mg/L)', pctsat='DO (% sat)', ddodt='dDO/dt (mg/L/d)')[y_var]))
  
  plot_out <- switch(
    style,
    'ggplot2' = {
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      
      . <- solar.time <- mod <- date <- col.mod <- col.obs <- obs <- '.ggplot.var'
      preds_ggplot <- v(DO_preds_all) %>%
        filter(as %in% y_var)
      if('conc' %in% names(y_lim)) {
        lim <- y_lim[['conc']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'conc' | (pure >= lim & mod >= lim & obs >= lim))
        lim <- y_lim[['conc']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'conc' | (pure <= lim & mod <= lim & obs <= lim))
      }
      if('pctsat' %in% names(y_lim)) {
        lim <- y_lim[['pctsat']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'pctsat' | (pure >= lim & mod >= lim & obs >= lim))
        lim <- y_lim[['pctsat']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'pctsat' | (pure <= lim & mod <= lim & obs <= lim))
      }
      if('ddodt' %in% names(y_lim)) {
        lim <- y_lim[['ddodt']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ddodt' | (pure >= lim & mod >= lim & obs >= lim))
        lim <- y_lim[['ddodt']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ddodt' | (pure <= lim & mod <= lim & obs <= lim))
      }
      g <- ggplot2::ggplot(preds_ggplot, ggplot2::aes(x=solar.time, group=date))
      # optional (only applies to sim models): 'pure' lines
      if(any(!is.na(preds_ggplot$pure))) g <- g + ggplot2::geom_line(ggplot2::aes(y=pure, color=col.pure), size=0.8, na.rm=TRUE)
      g + ggplot2::geom_point(ggplot2::aes(y=obs, color=col.obs), alpha=0.6, na.rm=TRUE) +
        ggplot2::geom_line(ggplot2::aes(y=mod, color=col.mod), size=0.8, na.rm=TRUE) +
        ggplot2::scale_color_identity(guide='none') +
        ggplot2::theme_bw() + 
        ggplot2::facet_grid(var ~ ., scales="free_y") + 
        ggplot2::xlab(params$xlab) + ggplot2::ylab(params$ylab)
    },
    'dygraphs' = {
      if(!requireNamespace("dygraphs", quietly=TRUE))
        stop("call install.packages('dygraphs') before plotting with style='dygraphs'")
      if(!requireNamespace("xts", quietly=TRUE))
        stop("call install.packages('xts') before plotting with style='dygraphs'")
      
      . <- '.dplyr.var'
      preds_xts <- v(DO_preds_all) %>%
        filter(as %in% y_var) %>%
        arrange(solar.time) %>%
        group_by(date) %>%
        do(., {
          out <- .[c(seq_len(nrow(.)),nrow(.)),]
          out[nrow(.)+1,c('pure','mod','obs')] <- NA
          out
        }) %>%
        ungroup()
      
      prep_dygraph <- function(y_var) { 
        prepped <- preds_xts %>% 
          filter(as==y_var) %>% 
          select(pure,mod,obs,solar.time) %>%
          mutate(solar.time=lubridate::force_tz(solar.time, Sys.getenv("TZ"))) %>% # dygraphs makes some funky tz assumptions. this seems to help.
          xts::xts(x=select(., -solar.time), order.by=.$solar.time, unique=FALSE, tzone=Sys.getenv("TZ"))
        if(all(is.na(prepped[,'pure']))) prepped <- prepped[,c('mod','obs')]
        prepped
      }
      if(length(y_var) > 1) {
        y_var <- y_var[1]
        warning("can only plot one dygraph y_var at a time for now; plotting ", y_var)
      }
      y_var_long <- preds_xts %>% filter(as==y_var) %>% slice(1) %>% .[['var']] %>% as.character()
      y_var_col <- params$colors[[y_var]]
      dat <- prep_dygraph(y_var)
      ymin <- max(c(min(c(unclass(dat)), na.rm=TRUE), y_lim[[y_var]][1]), na.rm=TRUE)
      ymax <- min(c(max(c(unclass(dat)), na.rm=TRUE), y_lim[[y_var]][2]), na.rm=TRUE)
      d <- dygraphs::dygraph(dat, xlab=params$xlab, ylab=y_var_long, group='plot_DO_preds')
      if(ncol(dat) == 3) d <- d %>% dygraphs::dySeries('pure', drawPoints = FALSE, label=paste0("Pure ", y_var_long), color=y_var_col[1])
      d %>% dygraphs::dySeries('mod', drawPoints = FALSE, label=paste0("Modeled ", y_var_long), color=y_var_col[2]) %>%
        dygraphs::dySeries('obs', drawPoints = TRUE, strokeWidth=0, label=paste0("Observed ", y_var_long), color=y_var_col[3]) %>%
        dygraphs::dyAxis('y', valueRange=(c(ymin,ymax)+(ymax-ymin)*c(-0.05,0.15))) %>%
        dygraphs::dyOptions(colorSaturation=1) %>%
        dygraphs::dyLegend(labelsSeparateLines = TRUE, width=300) %>%
        dygraphs::dyRangeSelector(height = 20)
    }
  )
  
  plot_out
}
