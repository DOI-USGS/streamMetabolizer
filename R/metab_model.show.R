#' @include metab_model-class.R
NULL

#' Display the metab_model object
#'
#' Print a metab_model object to the console.
#'
#' @param object metab_model to be displayed.
#' @importFrom utils head
#' @import dplyr
setMethod(
  "show", "metab_model",
  function(object) {
    cat("metab_model", if(class(object)[1] != "metab_model") paste0("of type ", class(object)[1]), "\n")
    if(!is.null(get_info(object))) {
      cat("  User-supplied metadata:\n")
      print(get_info(object))
    }
    cat("streamMetabolizer version", object@pkg_version, "\n")
    print(get_specs(object), header="Specifications:\n", prefix="  ")
    cat("Fitting time: ", tryCatch(get_fitting_time(object)[['elapsed']], error=function(e) "NA"), " secs elapsed\n", sep="")

    # set the sim.seed for a sec so the params and metab preds look the same.
    # setting the seed makes a copy of object, so we won't need to reset the
    # seed at the end of the function call
    seed <- get_specs(object)$sim_seed
    if(!is.null(seed) && is.na(seed)) object@specs$sim_seed <- Sys.time()

    # print parameter fits
    params <- get_params(object, uncertainty='ci', fixed='stars', messages=TRUE)
    if(!is.null(params)) {
      fixinfo <- if(any(grepl('\\*', params))) "(* = fixed value)" else ""
      cat("Parameters (", nrow(params), " date", if(nrow(params)!=1) "s", ")", fixinfo, ":\n", sep='')
      pretty_print_ddat(params, 'msgs.fit')
    }

    # print fitting warnings & errors
    warn_msgs <- summarize_stopwarn_msgs(params$warnings)
    stop_msgs <- summarize_stopwarn_msgs(params$errors)
    fit <- get_fit(object)
    if(!is.data.frame(fit) && is.list(fit)) {
      warn_msgs <- c(fit$warnings, warn_msgs)
      stop_msgs <- c(fit$errors, stop_msgs)
    }
    if(length(warn_msgs) > 0) cat("Fitting warnings:", paste0('  ', warn_msgs, collapse='\n'), sep='\n')
    if(length(stop_msgs) > 0) cat("Fitting errors:", paste0('  ', stop_msgs, collapse='\n'), sep='\n')

    # print metabolism predictions
    if(class(object)[1] == "metab_Kmodel") {
      return() # Kmodel doesn't do metabolism predictions
    }
    warn_msgs <- stop_msgs <- character(0)
    withCallingHandlers(
      tryCatch({
        metab_preds <- predict_metab(object)
        if(!is.null(metab_preds)) {
          cat("Predictions (", nrow(metab_preds), " date", if(nrow(metab_preds)!=1) "s", "):\n", sep='')
          pretty_print_ddat(metab_preds, 'msgs.pred')
          warn_msgs <- c(warn_msgs, summarize_stopwarn_msgs(metab_preds$warnings))
          stop_msgs <- c(stop_msgs, summarize_stopwarn_msgs(metab_preds$errors))
        }
      }, error=function(err) {
        stop_msgs <<- c(stop_msgs, err$message)
      }), warning=function(war) {
        warn_msgs <<- c(stop_msgs, war$message)
        invokeRestart("muffleWarning")
      })
    if(length(warn_msgs) > 0) cat("Prediction warnings:", paste0('  ', warn_msgs, collapse='\n'), sep='\n')
    if(length(stop_msgs) > 0) cat("Prediction errors:", paste0('  ', stop_msgs, collapse='\n'), sep='\n')
  }
)

#' Summarize a vector of warning or error messages
#'
#' Split ;-separated warning/error messages and condense into counts of each
#' unique message
#'
#' @keywords internal
summarize_stopwarn_msgs <- function(msgs) {
  if(is.null(msgs)) return(c())
  split_msgs <- unlist(strsplit(msgs, '; '))
  tbl_msgs <- sort(table(split_msgs))
  if(length(tbl_msgs) == 0) {
    c()
  } else {
    msgs_w_counts <- paste0(
      unname(tbl_msgs), " date", ifelse(unname(tbl_msgs)==1,"","s"), ": ",
      names(tbl_msgs))
    sort(msgs_w_counts)
  }
}

#' Compress warnings and errors in to a single column
#'
#' Compress two columns of warning and error messages into one short-hand column
#'
#' @param ddat a data.frame including warnings and errors columns
#' @param colname the name of the column where the summary should be placed
#' @param warnings.overall any general warnings (for the whole model) to be
#'   included in the summary
#' @param errors.overall any general errors (for the whole model) to be included
#'   in the summary
#' @import dplyr
#' @keywords internal
compress_msgs <- function(ddat, colname='messages', warnings.overall=c(), errors.overall=c()) {
  generr <- length(errors.overall) > 0 | grepl("overall errors", ddat$errors)
  genwarn <- length(warnings.overall) > 0 | grepl("overall warnings", ddat$warnings)
  dateerr <- gsub("overall errors(;)*", "", ddat$errors)
  datewarn <- gsub("overall warnings(;)*", "", ddat$warnings)
  messages <- errors <- warnings <- '.dplyr.var'
  ddat %>%
    mutate(
      messages = paste(
        # general errors and warnings
        ifelse(
          generr,
          ifelse(genwarn, 'w e', '  e'),
          ifelse(genwarn, 'w  ', '   ')),
        # date-specific errors and warnings
        ifelse(
          dateerr != '',
          ifelse(datewarn != '', 'W E', '  E'),
          ifelse(datewarn != '', 'W  ', '   '))
      )) %>%
    select(-warnings, -errors) %>%
    rename_with(~ colname, messages)
}

#' Format and print a summary of data frame of daily values
#'
#' Compress error and warning messages into one column and only print the first
#' few rows of data
#'
#' @import dplyr
#' @keywords internal
pretty_print_ddat <- function(ddat, msg.col) {
  if(!exists('warnings', ddat)) ddat$warnings <- NA
  if(!exists('errors', ddat)) ddat$errors <- NA
  ddat %>%
    head(10) %>%
    compress_msgs(colname=msg.col) %>%
    print()
  if(nrow(ddat) > 10) cat("  ...\n")
}
