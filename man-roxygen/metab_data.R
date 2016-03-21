#' @section Formatting \code{data}:
#' Unit-value model inputs passed via the \code{data} argument should
#' be formatted as a data.frame with column names and values that
#' depend on the model \code{type}, as follows.
#' (If all columns are optional, \code{data} may equal \code{NULL}.)
#' 
#' \describe{
#'   \item{\code{mle}, \code{bayes}, or \code{night}}{
#'     \tabular{lrrrrr}{
#'       \strong{solar.time} \tab \strong{DO.obs} \tab \strong{DO.sat} \tab \strong{depth} \tab \strong{temp.water} \tab \strong{light}\cr
#'       2050-03-14 15:10:00  \tab 10.1             \tab 14.2             \tab 0.5             \tab 21.8                 \tab 300.9          
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrr}{
#'       \code{solar.time         } \tab \code{DO.obs} \tab \code{DO.sat} \tab \code{depth} \tab \code{temp.water} \tab \code{light}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{10.1  } \tab \code{14.2  } \tab \code{0.5  } \tab \code{21.8      } \tab \code{300.9}
#'     }
#'   }
#'   \item{\code{Kmodel}}{
#'     \tabular{lrr}{
#'       \strong{solar.time} \tab \strong{discharge} \tab \strong{velocity}\cr
#'       2050-03-14 15:10:00  \tab 9                   \tab 2                 
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrr}{
#'       \code{solar.time         } \tab \code{discharge} \tab \code{velocity}\cr
#'       \code{2050-03-14 15:10:00} \tab \code{9        } \tab \code{2       }
#'     }
#'   }
#'   \item{\code{sim}}{
#'     \tabular{lrrrrr}{
#'       \strong{solar.time} \tab \strong{DO.obs} \tab \strong{DO.sat} \tab \strong{depth} \tab \strong{temp.water} \tab \strong{light}\cr
#'       2050-03-14 15:10:00  \tab 10.1             \tab 14.2             \tab 0.5             \tab 21.8                 \tab 300.9          
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
#'   \item{\code{bayes} or \code{night}}{
#'     \code{NULL}
#'   }
#'   \item{\code{mle}}{
#'     \tabular{lr}{
#'       \strong{date} \tab \strong{K600}\cr
#'       2050-03-14     \tab 5             
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lr}{
#'       \code{date      } \tab \code{K600}\cr
#'       \code{2050-03-14} \tab \code{5   }
#'     }
#'   }
#'   \item{\code{Kmodel}}{
#'     \tabular{lrrrrr}{
#'       \strong{date} \tab \strong{K600} \tab \strong{K600.lower} \tab \strong{K600.upper} \tab \strong{discharge.daily} \tab \strong{velocity.daily}\cr
#'       2050-03-14     \tab 5              \tab 4.5                  \tab 5.6                  \tab 9                         \tab 2                       
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrrr}{
#'       \code{date      } \tab \code{K600} \tab \code{K600.lower} \tab \code{K600.upper} \tab \code{discharge.daily} \tab \code{velocity.daily}\cr
#'       \code{2050-03-14} \tab \code{5   } \tab \code{4.5       } \tab \code{5.6       } \tab \code{9              } \tab \code{2             }
#'     }
#'   }
#'   \item{\code{sim}}{
#'     \tabular{lrrrr}{
#'       \strong{date} \tab \strong{DO.mod.1} \tab \strong{GPP} \tab \strong{ER} \tab \strong{K600}\cr
#'       2050-03-14     \tab 7.5                \tab 5             \tab 5            \tab 5             
#'     }
#'     
#'     \strong{Example}:
#'     \tabular{lrrrr}{
#'       \code{date      } \tab \code{DO.mod.1} \tab \code{GPP} \tab \code{ER} \tab \code{K600}\cr
#'       \code{2050-03-14} \tab \code{7.5     } \tab \code{5  } \tab \code{5 } \tab \code{5   }
#'     }
#'   }
#' }
