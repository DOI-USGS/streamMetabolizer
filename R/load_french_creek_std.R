#' Load a short dataset from French Creek using Bob Hall's code
#' 
#' This function requires \code{chron}, a package that is not formally required
#' by the \code{streamMetabolizer} package overall. Ensure that you have
#' \code{chron} installed or install it with \code{install.packages('chron')}.
#' 
#' This function produces a version of the French Creek dataset that agrees as 
#' much as possible with raw code from Bob Hall, for comparison to functions 
#' written in streamMetabolizer. Line numbers in comments (L22, etc.) refer to 
#' those in 'stream_metab_usa/starter_files/One station metab code.R'
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
    stop("chron package is needed for this function. Try install.packages('chron')")
  }
  french$dtime <- chron::chron(dates=as.character(french$date), times=as.character(french$time)) # L22. TZ is MST (L108, 142)
  tz_french <- lubridate::tz(convert_UTC_to_localtime(as.POSIXct("2012-09-10 00:00:00", tz="UTC"), latitude=41.33, longitude=-106.3, time.type="standard"))
  french$local.time <- with_tz(as.POSIXct(paste(french$date, french$time), format="%m/%d/%Y %H:%M:%S", tz="America/Denver"), tz_french) # need POSIXct for streamMetabolizer
  french$utc.time <- convert_localtime_to_UTC(french$local.time)
  french$solar.time <- convert_UTC_to_solartime(french$utc.time, longitude=-106.3, time.type='mean solar')
  
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
    # time <- time - 1/24 # if you adjust for DST here, also change longobs to (longobs-15) below. bob just goes with DST time for calculating jday, etc.
    jday<-as.numeric(trunc(time)-as.numeric(as.Date(year))) # bob does jday before longitude or DST adjustments. sM does after. should mostly affect midnight values, when it's dark anyway.
    E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
    LST<-as.numeric (time-trunc(time))
    ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440 # bob adjusts from local clock time. sM adjusts from UTC. since the earth rotates in 1436 mins rather than 1440, this causes a slight discrepancy (1-2 mins).
    solardel<- 23.439*sin(radi(360*((283+jday)/365)))
    hourangle<-(0.5-ST)*360
    theta<- acos(  sin(radi(solardel)) * sin(radi(lat)) +  cos(radi(solardel)) * cos(radi(lat)) *  cos(radi(hourangle)) )
    suncos<-ifelse(cos(theta)<0, 0, cos(theta))
    GI<- suncos*2326
    GI	
  }
  french$light <- lightest(time=french$dtime, lat=41.33,  longobs=106.3, longstd=90, year="2012-01-01") # L142

  # return w/ columns needed by streamMetabolzier
  oxy <- temp <- light <- ".dplyr.var"
  french <- dplyr::select(french, solar.time, DO.obs=oxy, DO.sat, depth, temp.water=temp, light=light)  
}

#' Generate outputs using Bob's code for comparison
#' 
#' Bob's code includes MLE and nighttime regression models. This function 
#' generates the output from those models, keeping the code as much intact as 
#' possible. The exception is that we're using solar.time rather than 
#' local.time, for consistency with streamMetabolizer's recommendations
#' 
#' This function requires the \code{chron} package, which is only suggested 
#' rather than required for the \code{streamMetabolizer} package. If you wish to
#' run this function, ensure that \code{chron} is installed or install it with 
#' \code{install.packages('chron')}.
#' 
#' @param french the French Creek dataset
#' @param K optional. If specified, a number for the K600 to assume (units of 
#'   1/d)
#' @param estimate character indicating the type of model to fit
#' @param start a character vector specifying the time at which the 'day' (the 
#'   time period to use in producing an estimate for a single date) begins. The 
#'   vector should have 2 elements, dates and times, to pass to chron()
#' @param end a character vector specifying the time at which the 'day' ends. 
#'   The vector should have 2 elements, dates and times, to pass to chron()
#' @param plot logical - should plots be produced?
#' @import dplyr
#' @importFrom unitted u v
#' @importFrom graphics abline plot points
load_french_creek_std_mle <- function(french, K=35, estimate=c('PRK','K','PR'), 
                                      start=c(dates="08/23/12", times="22:00:00"),
                                      end=c(dates="08/25/12", times="06:00:00"),
                                      plot=FALSE) {
  
  # require chron package
  if (!requireNamespace("chron", quietly = TRUE)) {
    stop("chron package is needed for this function. Try install.packages('chron')")
  }
  
  # define defaults for start and end here since they rely on the optional chron package
  if(missing(start)) 
    start <- chron::chron(dates="08/23/12", times="22:00:00") 
  else
    start <- do.call(chron::chron, as.list(start))
  if(missing(end)) 
    end <- chron::chron(dates="08/25/12", times="06:00:00") 
  else
    end <- do.call(chron::chron, as.list(end))
  
  # reformat, subset the data to just the requested day
  french <- v(french)
  french$dtime <- chron::chron(format(french$solar.time, "%m/%d/%y"), times=format(french$solar.time, "%H:%M:%S"))
  o2file <- french[french$dtime>=as.numeric(start) & french$dtime<=as.numeric(end), ]
  
  # set constants specific to French Creek data
  z=0.16
  bp=523
  ts=5/1440 # was 0.003422, which is ~= 5/1440=0.0034722
  
  # function to correct K600 to instantaneous KO2
  Kcor<-function (temp,K600) {
    #K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5 # bob & maite's original
    K600/(600/(1568 - 86.04*temp + 2.142*temp^2 - 0.0216*temp^3))^-0.5 # to match raymond et al. 2012
  }
  
  # function to calculate oxygen saturation given temperature and BP (mmHg). From Garcia and Gordon 1992 L&O
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
  }
  
  # fit the requested model, produce metabolism estimates
  estimate <- match.arg(estimate)
  ests <- switch(
    estimate,
    'PRK'={
      rivermetabK <- function(o2file, z, bp, ts){
        # likelihood function - returns the likelihood value for given GPP ER and K (which is vector MET)
        onestationmleK<-function(MET, temp, oxy, light, z, bp, ts) {
          # simulate DO
          metab<-numeric(length(temp))
          metab[1]<-oxy[1]
          for (i in 2:length(oxy)) {
            metab[i] <- metab[i-1]+
              ((MET[1]/z)*(light[i]/sum(light)))+ 
              MET[2]*ts/z+
              (Kcor(temp[i],MET[3]))*ts*(osat(temp[i],bp)-metab[i-1]) 
          }
          # calculate likelihood
          sqdiff<-(oxy-metab)^2 
          length(oxy)*(log(((sum(sqdiff)/length(oxy))^0.5)) +0.5*log(2*pi))   + ((2*sum(sqdiff)/length(oxy))^-1)*sum(sqdiff)
        }
        ##calculate metabolism by non linear minimization of MLE function
        river.mle<-nlm(onestationmleK, p=c(3,-5, 10), oxy=o2file$DO.obs, z=z, temp=o2file$temp.water, light=o2file$light, bp=bp, ts=ts)
        data.frame(GPP=river.mle$estimate[1], ER=river.mle$estimate[2], K=river.mle$estimate[3], lik=river.mle$minimum[1])
      }
      rivermetabK(o2file, z=z, bp=bp, ts=ts)
    },
    'K'={
      nightreg<-function(o2file, bp, ts){
        ##calculate delO/delt
        oxyf1<-stats::filter(o2file$DO.obs, rep(1/3,3), sides=2) # moving average on oxy data
        oxyf2<- oxyf1[c(-1,-length(oxyf1))] #trim the ends of the oxy data
        deltaO2<-((oxyf2[-1]-oxyf2[-length(oxyf2)])/ts)*1440
        
        #calc the do deficit
        temptrim <- o2file$temp.water[c(-2:-1,-nrow(o2file))] #Trim the first two and last one from the temp data to match the filter oxy data
        satdef<-osat(temptrim,bp)-oxyf2[-1]
        
        # calculate regression
        nreg<-lm(deltaO2~satdef)
        
        # plot
        if(plot) {
          plot(satdef,deltaO2)
          abline(nreg)
        }
        
        # extract and return coefficients + K600
        K600fromO2<-function (temp, KO2) { # function to estimate K600 from KO2  needed for nightime regression
          #((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2 # bob & maite's original
          ((600/(1568 - 86.04*temp + 2.142*temp^2 - 0.0216*temp^3))^-0.5)*KO2 # to match raymond et al. 2012
          
        }
        data.frame(GPP=as.numeric(NA), ER=NA, K=K600fromO2(mean(o2file$temp.water), coef(nreg)['satdef']), lik=as.numeric(NA))
      }
      nightreg(o2file, bp=bp, ts=ts*1440)
    },
    'PR'={
      rivermetab<-function(o2file, z, bp, ts, K){
        # likelihood function - returns the likelihood value for given GPP ER and K (which is vector MET)
        onestationmle<-function(MET,temp, oxy, light, z, bp, ts, K) {
          # simulate DO
          metab<-numeric(length(temp))
          metab[1]<-oxy[1]
          for (i in 2:length(oxy)) {
            metab[i]<-metab[i-1]+
              ((MET[1]/z)*(light[i]/sum(light)))+ 
              MET[2]*ts/z+
              (Kcor(temp[i],K))*ts*(osat(temp[i],bp)-metab[i-1]) 
          }
          ## calculate likelihood
          sqdiff<-(oxy-metab)^2 
          length(oxy)*(log(((sum(sqdiff)/length(oxy))^0.5)) +0.5*log(2*pi))   + ((2*sum(sqdiff)/length(oxy))^-1)*sum(sqdiff)
        }
        ##calculate metabolism by non linear minimization of MLE function
        river.mle<-nlm(onestationmle, p=c(3,-5), oxy=o2file$DO.obs, z=z, temp=o2file$temp.water, light=o2file$light, bp=bp, ts=ts, K=K)
        data.frame(GPP=river.mle$estimate[1], ER=river.mle$estimate[2], K=K, lik=river.mle$minimum[1])
      }
      rivermetab(o2file, z=z, bp=bp, ts=ts, K=K)
    })
  
  # plot data if PR or PRK (already plotted for K)
  if(plot) {
    if(estimate %in% c('PRK','PR')) {
      onestationplot<-function(GPP, ER, oxy, z, temp, K, light, bp, ts) {
        metab<-numeric(length(oxy))
        metab[1]<-oxy[1]
        for (i in 2:length(oxy)) { 
          metab[i] <- metab[i-1]+
            ((GPP/z)*(light[i]/sum(light)))+ 
            ER*ts/z+
            (Kcor(temp[i],K))*ts*(osat(temp[i],bp)-metab[i-1]) 
        }
        plot(seq(1:length(oxy)),metab, type="l",xlab="Time", ylab="Dissolved oxygen  (mg/L)", cex.lab=1.5, cex.axis=1.5, lwd=2 )
        points(seq(1:length(oxy)),oxy)
      }
      onestationplot(GPP=ests$GPP[1], ER=ests$ER[1], oxy=o2file$DO.obs, z=z, temp=o2file$temp.water, light=o2file$light, K=ests$K[1], bp=bp, ts=ts)
    }
  }
  
  # return estimates
  ests
}