#' Load a short dataset from French Creek using Bob Hall's code
#' 
#' Line numbers in comments (L22, etc.) refer to those in 
#' 'stream_metab_usa/starter_files/One station metab code.R'
#' 
#' @import dplyr
#' @importFrom unitted u
#' @importFrom utils read.csv
#' @importFrom lubridate with_tz
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_creek_std <- function(attach.units=TRUE) {
  # load the file
  file.name <- system.file("extdata", "french.csv", package="streamMetabolizer") # L10 data from French Creek, Hotchkiss and Hall, In press, Ecology
  french <- read.csv(file.name) 
  french <- french[french$station=="low",] # L12 subset to data from only one station (also, low has the cleanest data) (already subset in extdata)
  
  # remove NA oxys (1658) and remaining duplicates (n=1)
  french <- unique(french[!is.na(french$oxy),])
  
  # datetime
  if (!requireNamespace("chron", quietly = TRUE)) {
    stop("chron package is needed for this function. Try install.packages('chron')", call. = FALSE)
  }
  french$dtime <- chron::chron(dates=as.character(french$date), times=as.character(french$time)) # L22. TZ is MST (L108, 142)
  french$local.time <- with_tz(as.POSIXct(paste(french$date, french$time), format="%m/%d/%Y %H:%M:%S", tz="America/Denver"), "MST") # need POSIXct for streamMetabolizer
  
  # DO at sat
  osat <- function(temp, bp){ # L26-32
    ts<-log((298.15-temp) / (273.15 + temp))
    a0<-2.00907
    a1<-3.22014
    a2<-4.0501
    a3<-4.94457
    a4<- -0.256847
    a5<- 3.88767
    
    u<-10^(8.10765-(1750.286/(235+temp)))
    (exp(a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4+a5*ts^5))*((bp-u)/(760-u))*1.42905
  } # calculates oxygen saturation given temperature and BP (mmHg). From Garcia and Gordon 1992 L&O
  bpcalc <- function(bpst, alt, temp) { # L34-40
    bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+temp)))
  } # corrects for altitude.  From Colt's book. temp is degC, alt is m, and bpst is in inches of Hg. Temp is usually relative to a standard, 15 degC. 
  bpcalc(29.92,3150,15) # here's how to get very close to bp=523
  french$DO.sat <- osat(temp=french$temp, bp=523) # bp=523 mmHg for nightreg (L100), 595 for rivermetabK (L213) & rivermetab (L261). But 523 is probably correct (different from spring creek)
  
  # depth
  french$depth <- 0.16 # L261 gives 0.18, but Bob says it's actually 0.16
  
  # light
  lightest<- function (time, lat, longobs, longstd, year ) { # L107-135
    radi<-function(degrees){(degrees*pi/180)}
    time <- time - 1/24
    jday<-as.numeric(trunc(time)-as.numeric(as.Date(year)))
    E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
    LST<-as.numeric (time-trunc(time))
    ST<-LST+(3.989/1440)*(longstd-(longobs-15))+E/1440
    solardel<- 23.439*sin(radi(360*((283+jday)/365)))
    hourangle<-(0.5-ST)*360
    theta<- acos(  sin(radi(solardel)) * sin(radi(lat)) +  cos(radi(solardel)) * cos(radi(lat)) *  cos(radi(hourangle)) )
    suncos<-ifelse(cos(theta)<0, 0, cos(theta))
    GI<- suncos*2326
    GI	
  }
  french$light <- lightest(time=french$dtime, lat=41.33,  longobs=106.3, longstd=90, year="2012-01-01") # L142

  # return w/ columns needed by streamMetabolzier
  french <- dplyr::select(french, local.time, DO.obs=oxy, DO.sat, depth, temp.water=temp, light=light)  
}