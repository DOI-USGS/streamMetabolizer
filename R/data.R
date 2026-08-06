#' Example two-station (VFTS) input data
#'
#' A 30-day example dataset for fitting a two-station (upstream/downstream,
#' Variable Flow Two-Station) metabolism model with
#' \code{\link{metab_bayes_2s}}. It is a subset of the \code{VFTS-2}
#' (variable-travel-time) run from the published two-station metabolism modeling
#' dataset for a reach of the Colorado River in Glen Canyon, covering 2011-07-31
#' through 2011-08-29 at the source data's native 15-minute timestep, plus a
#' lead-in block of upstream DO observations (exactly as long as the longest
#' travel time in the dataset requires -- see \code{\link{metab_bayes_2s}}'s
#' "Two-station data requirements" section) immediately before 2011-07-31.
#'
#' @format A data.frame with 2904 rows and the 9 columns expected by
#'   \code{\link{metab_bayes_2s}}'s \code{data} argument, each carrying
#'   \code{\link[unitted]{unitted}} units matching \code{\link{mm_data}}:
#'   \describe{
#'     \item{solar.time}{POSIXct timestamp, UTC}
#'     \item{DO.obs.up}{dissolved oxygen observed at the upstream station,
#'       mgO2 L^-1}
#'     \item{DO.sat.up}{dissolved oxygen at equilibrium saturation at the
#'       upstream station, mgO2 L^-1}
#'     \item{DO.obs.down}{dissolved oxygen observed at the downstream
#'       station, mgO2 L^-1}
#'     \item{DO.sat.down}{dissolved oxygen at equilibrium saturation at the
#'       downstream station, mgO2 L^-1}
#'     \item{light}{within-day light proportion: light reaching the reach
#'       between each pair of upstream and downstream timepoints, summed over
#'       the upstream travel-time window and divided by that day's total
#'       light sum (Bishop et al. 2026 Eq. 2); unitless, range 0-1. Matches
#'       the quantity computed by \code{\link{mm_lag_light_2s}}, not the raw
#'       photosynthetically active radiation documented for \code{light} in
#'       \code{\link{mm_data}}.}
#'     \item{depth}{reach depth, m}
#'     \item{temp.water}{water temperature at the downstream station, degC}
#'     \item{travel.time}{reach travel time between the upstream and
#'       downstream stations, d}
#'   }
#'
#' @source Filtered to \code{model_run == 'VFTS-2'}; see
#'   \code{data-raw/two_station_example.R} for the extraction/renaming code.
#'
#'   Bishop, I.W., Deemer, B.R., Kennedy, T.A., Payn, R.A., Hall Jr, R.O. and
#'   Yackulic, C.B., 2026. A simplified two-station approach for modeling
#'   metabolism in dam tailwaters subject to diel flow variation. Limnology
#'   and Oceanography: Methods, p.e70066.
#'   \url{https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.1002/lom3.70066}
#'
#'   Data archived at ScienceBase:
#'   \url{https://www.sciencebase.gov/catalog/item/6887d457d4be024722b4aae2}
"two_station_example"

#' Raw, unformatted two-station (VFTS) input data
#'
#' The raw counterpart to \code{\link{two_station_example}}: upstream and
#' downstream sonde exports plus a modeled light series, but with no join, lag
#' alignment, or mean-solar-time conversion -- unlike
#' \code{\link{two_station_example}}, which has all of that processing
#' already done. \code{upstream}, \code{downstream}, and \code{light} are
#' NOT pre-aligned to a common time grid or to each other.
#'
#' @format A named list of three data.frames, each carrying
#'   \code{\link[unitted]{unitted}} units matching \code{\link{mm_data}}
#'   where a unit applies:
#'   \describe{
#'     \item{upstream}{86425 rows. \code{timestamp} (POSIXct, UTC, raw clock
#'       time -- not yet mean solar time); \code{DO.obs}, \code{DO.sat}
#'       (dissolved oxygen observed and at equilibrium saturation, mgO2
#'       L^-1)}
#'     \item{downstream}{90148 rows. \code{timestamp}; \code{DO.obs},
#'       \code{DO.sat} (mgO2 L^-1); \code{temp.water} (degC); \code{depth}
#'       (m, 165 NAs); \code{travel.time} (reach travel time, d)}
#'     \item{light}{101192 rows. \code{timestamp}; \code{light} (raw,
#'       per-timestep combined irradiance -- the source CSV's Direct and
#'       Indirect components summed, not the day-normalized,
#'       travel-time-weighted proportion that \code{\link{two_station_example}}'s
#'       \code{light} column holds; see \code{\link{mm_lag_light_2s}} for
#'       that later conversion. Confirmed as umol m^-2 s^-1 (PPFD).)}
#'   }
#'
#'   All three data.frames are truncated to the shared upstream/downstream
#'   deployment overlap, 2008-03-01 19:00:00 to 2011-01-19 20:45:00 UTC;
#'   downstream and light both continue for years past this window in their
#'   raw source files (downstream to 2014-02-28, light to 2025-01-01) but
#'   are cut here for consistency, since upstream data doesn't exist past
#'   2011-01-19.
#'
#' @source Upstream and downstream sonde data exported from the same
#'   Glen Canyon deployment described in \code{\link{two_station_example}};
#'   see \code{data-raw/two_station_raw_example.R} for the extraction/
#'   cleaning code.
"two_station_raw_example"
