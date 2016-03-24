#' Describe the requirements for an argument to metab()
#' 
#' @param type the type of model you want to fit
#' @param input the name of an argument to pass into metab()
#' @import dplyr
#' @importFrom unitted v u get_units
#' @examples 
#' metab_inputs('night','specs')
#' metab_inputs('bayes','data')
#' metab_inputs('Kmodel','data_daily')
#' metab_inputs('mle','info')
#' @export
metab_inputs <- function(type=c('bayes','mle','night','Kmodel','sim'), input=c('specs','data','data_daily','info')) {

  # check inputs
  type <- match.arg(type)
  input <- match.arg(input)
  
  if(input == 'specs') {
    paste0("specs(mm_name('",type,"'))", 
                 " # see ?mm_name, ?mm_specs for more options")
  } else if(input %in% c('data','data_daily')) {
    mfun <- paste0('metab_',type)
    eg <- eval(formals(mfun)[[input]])
    # reformat so there's a row each for units, format, example, and optional-T/F
    if(is.null(v(eg))) {
      'NULL'
    } else {
      . <- 'dplyr.var'
      data.frame(
        colname = {
          names(eg)
        }, 
        class = {
          sapply(unname(v(eg)), function(col) paste0(class(col), collapse=','))
        },
        units = {
          get_units(eg) %>%
            unname()
        },
        need = {
          opt <- attr(eg, 'optional')
          opt_vec <- if(opt[1]=='all') {
            rep('optional', length(eg))
          } else if(opt[1]=='none') {
            rep('required', length(eg))
          } else {
            ifelse(names(eg) %in% opt, 'optional', 'required') 
          }
        }
      )
    }
      # bind_rows %>%
      # u(c(type=NA, get_units(mm_data(everything())))[names(.)])
  } else if(input == 'info') {
    "info may be NULL, a list, or any other data you want to attach to the output of metab()"
  }
}