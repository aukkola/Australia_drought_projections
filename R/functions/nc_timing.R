read_nc_var <-  function(file, var) {
  nc <- nc_open(file)
  var <- ncvar_get(nc, var)
  nc_close(nc)
  return(var)
  
}


#This function is borrowed from "raster" package
#Needs to be this complicated because different models use different calendars
#and simple R functions can't deal with this
read_nc_time <- function(file) {
  dodays <- TRUE
  dohours <- FALSE
  doseconds <- FALSE
  
  nc <- nc_open(file)
  
  un <- ncatt_get(nc, "time")$unit	
  if (substr(un, 1, 10) == "days since") { 
    startDate = as.Date(substr(un, 12, 22))
  } else if (substr(un, 1, 11) == "hours since") { 
    dohours <- TRUE
    dodays <- FALSE
    startTime <- substr(un, 13, 30)
    mult <- 3600
  } else if (substr(un, 1, 13) == "seconds since") { 
    doseconds <- TRUE
    dodays <- FALSE
    startTime = as.Date(substr(un, 15, 31))
    mult <- 1
  } else if (substr(un, 1, 12) == "seconds from") { 
    doseconds <- TRUE
    dodays <- FALSE
    startTime = as.Date(substr(un, 14, 31))
    mult <- 1
  } else {
    return(x)
  }
  if (!dodays) {
    start <- strptime(startTime, "%Y-%m-%d %H:%M:%OS", tz = "UTC")
    if (is.na(start)) start <- strptime(startTime, "%Y-%m-%d", tz = "UTC")
    if (is.na(start)) return(x)
    startTime <- start
    time <- startTime + as.numeric(getZ(x)) * mult
    time <- as.character(time)
    if (!is.na(time[1])) {
      x@z <- list(time)
      names(x@z) <- as.character('Date/time')
    }	
  } else if (dodays) {
    # cal = nc$var[[zvar]]$dim[[dim3]]$calendar ?
    cal <- ncdf4::ncatt_get(nc, "time", "calendar")
    if (! cal$hasatt ) {
      greg <- TRUE
    } else {
      cal <- cal$value
      if (cal =='gregorian' | cal =='proleptic_gregorian' | cal=='standard') {
        greg <- TRUE
      } else if (cal == 'noleap' | cal == '365 day' | cal == '365_day') { 
        greg <- FALSE
        nday <- 365
      } else if (cal == '360_day') { 
        greg <- FALSE
        nday <- 360
      } else {
        greg <- TRUE
        warning('assuming a standard calender:', cal)
      }
    }
    time <- ncvar_get(nc, "time")
    if (greg) {
      time <- as.Date(time, origin=startDate)
    } else {
      startyear <-  as.numeric( format(startDate, "%Y") )
      startmonth <- as.numeric( format(startDate, "%m") )
      startday <- as.numeric( format(startDate, "%d") )
      year <- trunc( as.numeric(time)/nday )
      doy <- (time - (year * nday))
      origin <- paste(year+startyear, "-", startmonth, "-", startday, sep='')
      time <- as.Date(doy, origin=origin)		
    }
 
  }
  nc_close(nc)
  
  return(time)
}

