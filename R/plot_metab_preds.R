#' Plot predictions produced with predict_DO
#' 
#' Plots modeled values as lines, observed values as points
#' 
#' @param metab_preds a data.frame of predictions such as that returned by 
#'   predict_metab()
#' @param y_var character. Should the plot display predicted values of GPP, ER, 
#'   and/or K600? The default is to plot all three.
#' @param style character indicating which graphics package to use
#' @param y_lim list of named vectors, each of which has length 2 and is numeric
#'   and has a name in the possible values of y_var. NA within a vector
#'   indicates that the data range should be used. for ggplot2, y_lim is only
#'   used to exclude values outside that range and is ignored if the data span a
#'   narrower range
#' @examples 
#' \dontrun{
#' mm <- metab_night(v(french))
#' plot_metab_preds(predict_DO(mm)[1:360,])
#' }
#' @import dplyr
#' @importFrom unitted v
#' @export
plot_metab_preds <- function(metab_preds, y_var=c('GPP','ER','K600'), 
                          style=c('ggplot2','dygraphs'),
                          y_lim=list(GPP=c(NA,NA), ER=c(NA,NA), K600=c(NA,NA))) {
 
  style <- match.arg(style)
  y_var <- match.arg(y_var, several.ok=TRUE)
  
  params <- list(
    xlab='Local date',
    ylab='Predictions',
    colors=list(GPP=c('#A64B00','#FF7400'), ER=c('#007929','#23BC47'), K600=c('#05326D','#4282D3'))
  )
  
  metab.mod <- '.dplyr.var'
  metab_preds_GPP <- mutate(
    metab_preds, as='GPP', val=GPP, var='GPP (g m^-2 d^-1)', col1=params$colors[['GPP']][1], col2=params$colors[['GPP']][2], lab='GPP~(g~m^-2~d^-1)')
  metab_preds_ER <- mutate(
    metab_preds, as='ER', val=ER, var='ER (g m^-2 d^-1)', col1=params$colors[['ER']][1], col2=params$colors[['ER']][2], lab='ER~(g~m^-2~d^-1)')
  metab_preds_K600 <- mutate(
      metab_preds, as='K600', val=K600, var='K600 (d^-1)', col1=params$colors[['K600']][1], col2=params$colors[['K600']][2], lab='K600~(d^-1)')
  
  metab_preds_all <- bind_rows(metab_preds_GPP, metab_preds_ER, metab_preds_K600) %>%
    mutate(var=ordered(var, c(GPP='GPP (g m^-2 d^-1)', ER='ER (g m^-2 d^-1)', K600='K600 (d^-1)')[y_var]))
  
  switch(
    style,
    'ggplot2' = {
      if(!requireNamespace("ggplot2", quietly=TRUE))
        stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
      
      . <- local.date <- col1 <- col2 <- '.ggplot.var'
      preds_ggplot <- v(metab_preds_all) %>%
        filter(as %in% y_var)
      if('GPP' %in% names(y_lim)) {
        lim <- y_lim[['GPP']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'GPP' | (val >= lim))
        lim <- y_lim[['GPP']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'GPP' | (val <= lim))
      }
      if('ER' %in% names(y_lim)) {
        lim <- y_lim[['ER']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ER' | (val >= lim))
        lim <- y_lim[['ER']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'ER' | (val <= lim))
      }
      if('K600' %in% names(y_lim)) {
        lim <- y_lim[['K600']][1]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'K600' | (val >= lim))
        lim <- y_lim[['K600']][2]; if(!is.na(lim)) preds_ggplot <- filter(preds_ggplot, as != 'K600' | (val <= lim))
      }
      g <- ggplot2::ggplot(preds_ggplot, ggplot2::aes(x=local.date)) +
        ggplot2::geom_line(ggplot2::aes(y=val, color=col1)) +
        ggplot2::geom_point(ggplot2::aes(y=val, color=col1)) +
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
      
      stop("no dygraphs option yet")
    }
  )
}