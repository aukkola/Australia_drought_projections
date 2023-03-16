library(raster)
library(zoo)
library(ncdf4)
library(lubridate)

#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"

#Preindustrial baselines
baseline_start <- c(1910)#, 1951)
baseline_end   <- c(1945)#, 1910) 


#End and start years (only save data for 1902-2100 to match the period common
#to all models and obs)
time_period <- c(1910, 2019)


#Global warming levels
envelopes <- c(1,2,3)


#Output directory
outdir <- paste0(path, "/Global_warming_levels/baseline_", baseline_start,
                 "_", baseline_end, "/observed/")
dir.create(outdir)


#Observed temp file
obs_file <- paste0(path, "/Observed_tair/Obs_mean_time_series/", 
                   "Monthly_mean_temperature_time_series_Aus.nc")
                   
nc <- nc_open(obs_file)

#variable called tmax because of a cdo command, but the file is tmean
#couldn't be bothered to fix variable name
data <- ncvar_get(nc, "tmax")


time <- ncvar_get(nc, "time")

time_units <- strsplit(ncatt_get(nc, "time")$units, "days since ")[[1]][2]

#Convert to Y-M-D h-m-s
time_date <- as.POSIXct(time*86400, origin=time_units, tz="GMT")

#Get years
years <- year(time_date)


#select years
start_ind <- which(years == time_period[1])[1]
end_ind   <- tail(which(years == time_period[2]), n=1)


monthly_temp <- data[start_ind:end_ind]


#Calculate annual temp




#Calculate 10yr running mean 
mean_10yr <- rollmean(annual_temp, k=10, fill=NA)

#Dataset years
data_yrs = data$Time[start_ind:end_ind]
      
for (b in 1:length(baseline_start)) {
  
  #Calculate pre-industrial mean for baseline
  preind_mean <- mean(mean_10yr[which(data_yrs >= baseline_start[b] & data_yrs <= baseline_end[b])],
                      na.rm=TRUE)
  
  
  yr_ind <- lapply(envelopes, function(x) which(mean_10yr >= (preind_mean+x-0.3) &
                                                  mean_10yr <=  (preind_mean+x+0.3)))
      
  
  #Convert data to monthly indices
  
  #Start year for decadal slice (yr_ind minus 9 years)
  start_ind <- lapply(yr_ind, function(x) x-9)
  
  #Then convert to monthly indices
  #awful mess but takes start and end index for each decadal slice and converts to monthly
  monthly_indices <- lapply(1:length(start_ind), function(x) if (length(start_ind[[x]]) > 0)
    ((min(start_ind[[x]])-1)*12+1) :  #previous year plus one month
      (max(yr_ind[[x]])*12)) #end year times 12 months
  
  names(monthly_indices) <- paste0(envelopes, "deg")
  
}


outfile <- paste0(outdir, "/Monthly_indices_global_warming_levels_", min(envelopes), 
                  "-", max(envelopes), "deg_CRUTEM5.0_observed.rds")

#Save as R object
saveRDS(monthly_indices, outfile)
        


