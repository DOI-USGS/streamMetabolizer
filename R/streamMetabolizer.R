#' Functions for calculating ecosystem metabolism in streams
#' 
#' This package contains functions that help the user prepare data and model
#' metabolism in streams and rivers.
#' 
#' @section Creating new columns:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{calc_depth}}
#'   
#'   \item \code{\link{calc_DO_at_sat}}
#'   
#'   \item \code{\link{calc_DO_deficit}}
#'   
#'   \item \code{\link{calc_is_daytime}}
#'   
#'   \item \code{\link{calc_solar_insolation}}
#'   
#'   \item \code{\link{calc_sun_rise_set}}
#'   
#'   }
#'   
#' @section Converting existing data:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{convert_date_to_doyhr}}
#'   
#'   \item \code{\link{convert_localtime_to_UTC}}
#'   
#'   \item \code{\link{convert_UTC_to_solartime}}
#'   
#'   \item \code{\link{convert_k600_to_kGAS}}
#'   
#'   \item \code{\link{convert_PAR_to_SW}}
#'   
#'   }
#'   
#' @section Model metabolism:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{mm_name}} 1. Choose a model structure
#'   
#'   \item \code{\link{specs}} 2. Set the specifications
#'   
#'   \item \code{\link{metab}} 3. Fit the model
#'   
#'   }
#'   
#' @section Inspect model results:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{predict_metab}}
#'   
#'   \item \code{\link{predict_DO}}
#'   
#'   \item \code{\link{plot_metab_preds}}
#'   
#'   \item \code{\link{plot_DO_preds}}
#'   
#'   \item \code{\link{get_fit}}
#'   
#'   \item \code{\link{get_specs}}
#'   
#'   \item \code{\link{get_data}}
#'   
#'   \item \code{\link{get_data_daily}}
#'   
#'   \item \code{\link{get_info}}
#'   
#'   \item \code{\link{get_mcmc}} (Bayesian models only)
#'   
#'   \item \code{\link{get_fitting_time}}
#'   
#'   \item \code{\link{get_version}}
#'   
#'   }
#'   
#' @docType package
#' @name streamMetabolizer
#' @aliases streamMetabolizer-package
NULL