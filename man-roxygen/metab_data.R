#' @section Formatting \code{data}:
#' Unit-value model inputs passed via the \code{data} argument should
#' be formatted as a data.frame with column names and values that
#' depend on the model \code{type}, as follows.
#' (If all columns are optional, \code{data} may equal \code{NULL}.)
#' 
#' \describe{
#'   \item{\code{mle} or \code{night}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       solar.time        \tab POSIXct,POSIXt  \tab                 \tab required      \cr
#'       DO.obs            \tab numeric         \tab mgO2 L^-1       \tab required      \cr
#'       DO.sat            \tab numeric         \tab mgO2 L^-1       \tab required      \cr
#'       depth             \tab numeric         \tab m               \tab required      \cr
#'       temp.water        \tab numeric         \tab degC            \tab required      \cr
#'       light             \tab numeric         \tab umol m^-2 s^-1  \tab required      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrr}{
#'       \code{solar.time         } \tab \code{DO.obs} \tab \code{DO.sat} \tab \code{depth} \tab \code{temp.water} \tab \code{light}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{10.1  } \tab \code{14.2  } \tab \code{0.5  } \tab \code{21.8      } \tab \code{300.9}
#'     }
#'   }
#'   \item{\code{bayes}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       solar.time        \tab POSIXct,POSIXt  \tab                 \tab required      \cr
#'       DO.obs            \tab numeric         \tab mgO2 L^-1       \tab required      \cr
#'       DO.sat            \tab numeric         \tab mgO2 L^-1       \tab required      \cr
#'       depth             \tab numeric         \tab m               \tab required      \cr
#'       temp.water        \tab numeric         \tab degC            \tab required      \cr
#'       light             \tab numeric         \tab umol m^-2 s^-1  \tab required      \cr
#'       discharge         \tab numeric         \tab m^3 s^-1        \tab optional      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrrr}{
#'       \code{solar.time         } \tab \code{DO.obs} \tab \code{DO.sat} \tab \code{depth} \tab \code{temp.water} \tab \code{light} \tab \code{discharge}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{10.1  } \tab \code{14.2  } \tab \code{0.5  } \tab \code{21.8      } \tab \code{300.9} \tab \code{9        }
#'     }
#'   }
#'   \item{\code{Kmodel}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       solar.time        \tab POSIXct,POSIXt  \tab                 \tab optional      \cr
#'       discharge         \tab numeric         \tab m^3 s^-1        \tab optional      \cr
#'       velocity          \tab numeric         \tab m s^-1          \tab optional      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrr}{
#'       \code{solar.time         } \tab \code{discharge} \tab \code{velocity}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{9        } \tab \code{2       }
#'     }
#'   }
#'   \item{\code{sim}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       solar.time        \tab POSIXct,POSIXt  \tab                 \tab required      \cr
#'       DO.obs            \tab numeric         \tab mgO2 L^-1       \tab optional      \cr
#'       DO.sat            \tab numeric         \tab mgO2 L^-1       \tab required      \cr
#'       depth             \tab numeric         \tab m               \tab required      \cr
#'       temp.water        \tab numeric         \tab degC            \tab required      \cr
#'       light             \tab numeric         \tab umol m^-2 s^-1  \tab required      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrr}{
#'       \code{solar.time         } \tab \code{DO.obs} \tab \code{DO.sat} \tab \code{depth} \tab \code{temp.water} \tab \code{light}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{10.1  } \tab \code{14.2  } \tab \code{0.5  } \tab \code{21.8      } \tab \code{300.9}
#'     }
#'   }
#' }
#' @section Formatting \code{data_daily}:
#' Daily-value model inputs passed via the \code{data_daily} argument should
#' be formatted as a data.frame with column names and values that
#' depend on the model \code{type}, as follows.
#' (If all columns are optional, \code{data_daily} may equal \code{NULL}.)
#' 
#' \describe{
#'   \item{\code{night}}{
#'     \code{NULL}
#'   }
#'   \item{\code{mle}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       date              \tab Date            \tab                 \tab optional      \cr
#'       K600              \tab numeric         \tab d^-1            \tab optional      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lr}{
#'       \code{date      } \tab \code{K600}\cr
#'       \code{2050-03-14} \tab \code{5   }
#'     }
#'   }
#'   \item{\code{bayes}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       date              \tab Date            \tab                 \tab optional      \cr
#'       discharge.daily   \tab numeric         \tab m^3 s^-1        \tab optional      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lr}{
#'       \code{date      } \tab \code{discharge.daily}\cr
#'       \code{2050-03-14} \tab \code{9              }
#'     }
#'   }
#'   \item{\code{Kmodel}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       date              \tab Date            \tab                 \tab required      \cr
#'       K600              \tab numeric         \tab d^-1            \tab required      \cr
#'       K600.lower        \tab numeric         \tab d^-1            \tab optional      \cr
#'       K600.upper        \tab numeric         \tab d^-1            \tab optional      \cr
#'       discharge.daily   \tab numeric         \tab m^3 s^-1        \tab optional      \cr
#'       velocity.daily    \tab numeric         \tab m s^-1          \tab optional      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrr}{
#'       \code{date      } \tab \code{K600} \tab \code{K600.lower} \tab \code{K600.upper} \tab \code{discharge.daily} \tab \code{velocity.daily}\cr
#'       \code{2050-03-14} \tab \code{5   } \tab \code{4.5       } \tab \code{5.6       } \tab \code{9              } \tab \code{2             }
#'     }
#'   }
#'   \item{\code{sim}}{
#'     \tabular{llll}{
#'       \strong{colname} \tab \strong{class} \tab \strong{units} \tab \strong{need}\cr
#'       date              \tab Date            \tab                 \tab required      \cr
#'       DO.mod.1          \tab numeric         \tab mgO2 L^-1       \tab optional      \cr
#'       GPP               \tab numeric         \tab gO2 d^-1 m^-2   \tab required      \cr
#'       ER                \tab numeric         \tab gO2 d^-1 m^-2   \tab required      \cr
#'       K600              \tab numeric         \tab d^-1            \tab required      
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrr}{
#'       \code{date      } \tab \code{DO.mod.1} \tab \code{GPP} \tab \code{ER} \tab \code{K600}\cr
#'       \code{2050-03-14} \tab \code{7.5     } \tab \code{5  } \tab \code{5 } \tab \code{5   }
#'     }
#'   }
#' }
