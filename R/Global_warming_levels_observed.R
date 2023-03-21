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
time_period <- c(1910, 2022)


#Global warming levels
envelopes <- c(1,2,3)


#Output directory
outdir <- paste0(path, "/Global_warming_levels/baseline_", baseline_start,
                 "_", baseline_end, "/observed/")
dir.create(outdir, recursive=TRUE)


#Observed temp file
obs_file <- paste0(path, "/Bom_obs_temp_data/BoM_Aus_mean_temperature_", 
                   "time_series_downloaded_17Mar2023.txt")
                   
#Read data
data <- read.table(obs_file)     
                   
#Geat years
data_yrs <- as.numeric( substring(data[,1], 1,4 ))


#select years
start_ind <- which(data_yrs == time_period[1])
end_ind   <- which(data_yrs == time_period[2])



#Calculate 10yr running mean 
mean_10yr <- rollmean(data[,2], k=10, fill=NA)



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
                  "-", max(envelopes), "deg_BoM_observed.rds")

#Save as R object
saveRDS(monthly_indices, outfile)
        


