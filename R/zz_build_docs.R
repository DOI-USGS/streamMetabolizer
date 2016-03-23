#' @include metab_inputs.R metab_bayes.R metab_mle.R metab_night.R metab_sim.R metab_Kmodel.R mm_data.R
NULL

#' Format a data.frame for inclusion in a roxygen header
#' 
#' Modified from Hadley Wickham's function at http://r-pkgs.had.co.nz/man.html
#' 
#' @keywords internal
zz_tabular <- function(df, bold_headers=TRUE, code=FALSE, ...) {
  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- vapply(df, align, character(1))
  
  cols <- 
    mapply(
      function(colname, colvec) {
        c(if(bold_headers) paste0("\\strong{", colname, "}") else colname,
          #paste(rep('-', nchar(colname)), collapse=''), 
          as.character(colvec)) 
      }, colname=as.list(names(df)), colvec=df) %>%
    as.data.frame() %>%
    lapply(format, ...)
  
  if(code) {
    cols <- lapply(cols, function(col) 
      paste0("\\code{", col, "}"))
  }
  
  cols <- as_data_frame(cols)
  
  contents <- do.call(
    "paste", 
    c(cols, list(sep = " \\tab ", collapse = "\\cr\n  ")))
  
  . <- 'dplyr.var'
  paste(
    "\\tabular{", 
    paste(col_align, collapse = ""), 
    "}{\n  ", 
    contents, 
    "\n}\n", 
    sep = "") %>%
    strsplit('\n') %>%
    .[[1]]
}

#' Generate doc text for the `metab()` documentation
#' 
#' Results get written to man-roxygen/metab_data.R
#' 
#' @keywords internal
zz_build_docs <- function() {
  
  . <- 'dplyr.var'
  c("@section Formatting \\code{data}:",
    "Unit-value model inputs passed via the \\code{data} argument should",
    "be formatted as a data.frame with column names and values that",
    "depend on the model \\code{type}, as follows.",
    "(If all columns are optional, \\code{data} may equal \\code{NULL}.)",
    "",
    "\\describe{",
    c(paste0("  \\item{\\code{mle} or \\code{night}}{"),
      paste0("    ", c(
        zz_tabular(metab_inputs('mle', 'data')),
        "",
        "\\strong{Example}:",
        zz_tabular(v(eval(formals(metab_mle)$data)), bold_headers=FALSE, code=TRUE)
      )),
      "  }"),
    do.call(c, lapply(c('bayes','Kmodel','sim'), function(type) {
      c(paste0("  \\item{\\code{",type,"}}{"),
        paste0("    ", c(
          zz_tabular(metab_inputs(type, 'data')),
          "",
          "\\strong{Example}:",
          zz_tabular(v(eval(formals(paste0("metab_",type))$data)), bold_headers=FALSE, code=TRUE)
        )),
        "  }")
    })),
    "}",
    
    "@section Formatting \\code{data_daily}:",
    "Daily-value model inputs passed via the \\code{data_daily} argument should",
    "be formatted as a data.frame with column names and values that",
    "depend on the model \\code{type}, as follows.",
    "(If all columns are optional, \\code{data_daily} may equal \\code{NULL}.)",
    "",
    "\\describe{",
    c(paste0("  \\item{\\code{night}}{"),
      paste0("    ", "\\code{", 
             metab_inputs('night', 'data_daily'),
             "}"),
      "  }"),
    do.call(c, lapply(c('mle','bayes','Kmodel','sim'), function(type) {
      c(paste0("  \\item{\\code{",type,"}}{"),
        paste0("    ", c(
          zz_tabular(metab_inputs(type, 'data_daily')),
          "",
          "\\strong{Example}:",
          zz_tabular(v(eval(formals(paste0("metab_",type))$data_daily)), bold_headers=FALSE, code=TRUE)
        )),
        "  }")
    })),
    "}"
  ) %>% 
    paste0("#' ", .) %>%
    writeLines('man-roxygen/metab_data.R')
  
  invisible()
}
zz_build_docs()