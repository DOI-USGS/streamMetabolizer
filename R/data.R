#' Example two-station (VFTS) input data
#'
#' A six-year example dataset for fitting a two-station (upstream/downstream,
#' Variable Flow Two-Station) metabolism model with
#' \code{\link{metab_bayes_2s}}. It is the complete \code{VFTS-2}
#' (variable-travel-time) run from the published two-station metabolism modeling
#' dataset for a reach of the Colorado River in Glen Canyon, covering
#' 2008-03-11 through 2014-02-28 at the source data's native 15-minute
#' timestep.
#'
#' The record is shipped whole rather than trimmed to a clean stretch, so it
#' still contains the 90 interruptions the sensors actually produced -- every
#' one of them a full day or longer. Combined with the 06:00-06:00 day window
#' (two-station days do not run midnight to midnight), this means the 1748
#' calendar dates spanned here yield 1551 modeled days: 91 dates do not fill a
#' complete 24-hour window, and a further 15 hold gaps in their upstream
#' oxygen columns. Days dropped for either reason are reported as
#' \code{valid_day=FALSE} rather than silently omitted. That messiness is the
#' point -- it exercises the day-validity handling that real sensor records
#' demand, which a gap-free excerpt would leave untested.
#'
#' The upstream and downstream oxygen columns are contemporaneous: the
#' \code{DO.obs.up} value on a given row was measured at that row's own
#' \code{solar.time}, not at the earlier time the water left the upstream
#' station. Pairing each downstream observation with the upstream water that
#' became it is done during fitting, using \code{travel.time}; each modeled
#' row therefore reaches back to an earlier row for its upstream value, and
#' rows with no such row available (at the very start of the record, and
#' immediately after each gap) are not modeled. Stored data must not have that
#' shift already applied.
#'
#' @format A data.frame with 159072 rows and the 9 columns expected by
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
#'   Fitting all 1551 days at once is a long computation; subset to a handful
#'   of days for interactive use.
#'
#' @source Drawn from the \code{VFTS-2} (variable-travel-time) run, except for
#'   the two upstream oxygen columns, which come from the same file's
#'   \code{VFTS-3} run. Both runs cover the same reach and the same physical
#'   sensors, but \code{VFTS-3} was run at a single fixed travel time, and only
#'   its constant shift can be undone exactly to recover the contemporaneous
#'   upstream series this dataset needs. See
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
#' downstream sonde exports plus a modeled light series, but with no join,
#' no gap filling, no light conversion, and no mean-solar-time conversion --
#' unlike \code{\link{two_station_example}}, which has all of that processing
#' already done. \code{upstream}, \code{downstream}, and \code{light} are
#' NOT pre-aligned to a common time grid or to each other.
#'
#' Neither dataset has the upstream travel-time shift applied: in both, the
#' upstream and downstream oxygen values at a given timestamp are
#' contemporaneous observations. That shift belongs to model fitting, which
#' applies it from \code{travel.time}, so applying it to stored data would
#' double it.
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
