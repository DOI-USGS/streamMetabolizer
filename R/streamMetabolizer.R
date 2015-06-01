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
#'   \item \code{\link{calc_schmidt}}
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
#'   \item \code{\link{convert_GMT_to_solartime}}
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
#'   \item \code{\link{metab_mle}}
#'   
#'   }
#'   
#' @section Internal functions:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{metab_model-class}}
#'   
#'   \item \code{\link{metab_model_interface}}
#'   
#'   \item \code{\link{mm_data}}
#'   
#'   \item \code{\link{mm_validate_data}}
#'   
#'   }
#'   
#' @docType package
#' @name streamMetabolizer
#' @aliases streamMetabolizer-package
NULL