library(raster)
library(zoo)
library(ncdf4)

#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"

source(paste0(path, "/scripts/R/functions/nc_timing.R"))


#Experiments
experiments <- list.files(paste0(path, "/Global_mean_temp/CABLE_AWRA/"), pattern="rcp")

#Only save indices up to 2100 (some models go well beyond this but want to keep
#this consistent)
end_year <- 2099


#Preindustrial baselines (to match CABLE/AWRA data availability)
baseline_start <- c(1970)
baseline_end   <- c(2005) 


#Global warming levels
#these corresponds to 2 and 3 degree worlds really as using a late baseline by 
#which about 1deg of warming has already occurred
envelopes <- c(1,2) 

#Output directory
outdir <- paste0(path, "/Global_warming_levels/CABLE_AWRA/")


#Loop through experiments
for (exp in 1:length(experiments)) {
  
  #Find bias correction methods
  bc_methods <- list.files(paste0(path, "/Global_mean_temp/CABLE_AWRA/", experiments[exp]))
  
  for (b in 1:length(bc_methods)) {
    
    #List models
    models <- list.files(paste0(path, "/Global_mean_temp/CABLE_AWRA/", experiments[exp], 
                                "/", bc_methods[b]))
    
    #Loop through models
    for (m in 1:length(models)) {
      
    
      #Output file
      outfile <- paste0("Monthly_indices_global_warming_levels_", min(envelopes), 
                        "-", max(envelopes), "deg_", experiments[exp], "_", bc_methods[b],
                        "_", models[m], ".rds")
      
      if (file.exists(outfile)) next
      
      
      #Get historical data 
      hist_file <- list.files(paste0(path, "/Global_mean_temp/CABLE_AWRA/historical/", 
                                     bc_methods[b], "/", models[m]), full.names=TRUE)
      
      hist_data <- read_nc_var(hist_file, var="tasmin") #called tasmin because of CDO combining tasmin and tasmax

      
      
      #Get future data
      fut_file <- list.files(paste0(path, "/Global_mean_temp/CABLE_AWRA/",
                             experiments[exp], "/",  bc_methods[b], "/", models[m]),
                             full.names=TRUE)
      
      fut_data <- read_nc_var(fut_file, var="tasmin") #called tasmin because of CDO combining tasmin and tasmax
      
      
      #Combine historical and future periods
      all_data <- c(hist_data, fut_data)
      
      
      
      #Calculate mean annual temperature
      ind <- seq(1, by=12, length.out=length(all_data)/12)
      
      annual_temp <- sapply(ind, function(x) mean(all_data[x:(x+11)]))
      
      
      #Get data years
      start_yr <- format(read_nc_time(hist_file), format="%Y")[1]
      
      
      if (start_yr != 1960) {
        stop("wrong historical start year")
      }
      
      data_yrs <- seq(start_yr, by=1, length.out=length(annual_temp))

      
      #Extract time period to stop in 2100
      #(some models end in 2099, use this instead)
      if (max(data_yrs) < end_year) {
        end <- length(data_yrs)
      } else {
        end   <- which(data_yrs == end_year)
      }

      data_yrs    <- data_yrs[1:end]
      annual_temp <- annual_temp[1:end]


      
      #Calculate 10yr running mean 
      mean_10yr <- rollmean(annual_temp, k=10, fill=NA)
      
      
      
      for (b in 1:length(baseline_start)) {
        
        #Calculate pre-industrial mean for baseline
        preind_mean <- mean(mean_10yr[which(data_yrs >= baseline_start[b] & data_yrs <= baseline_end[b])],
                            na.rm=TRUE)
        
        
        yr_ind <- lapply(envelopes, function(x) which(mean_10yr >= (preind_mean+x-0.3) &
                                                      mean_10yr <=  (preind_mean+x+0.3)))

        
        
        #Take out the historical period from the indices so they are relative to
        #the start of the future period
        
        #yr_ind_rebased <- lapply(yr_ind, function(x) x- (which(data_yrs == 2015)-1))
        
        
        
        # 
        # #Check that none of the time slices are during the historical period
        # if (any (unlist(yr_ind_rebased) > 1)) {
        #   stop("index during historical period")
        # }
        # 
        
        
        #Convert data to monthly indices
        
        #Start year for decadal slice (yr_ind minus 9 years)
        start_ind <- lapply(yr_ind, function(x) x-9)
        
        #Then convert to monthly indices
        #awful mess but takes start and end index for each decadal slice and converts to monthly
        monthly_indices <- lapply(1:length(start_ind), function(x) if (length(start_ind[[x]]) > 0)
                                                                    ((min(start_ind[[x]])-1)*12+1) :  #previous year plus one month
                                                                    (max(yr_ind[[x]])*12)) #end year times 12 months
                                    
        names(monthly_indices) <- paste0(envelopes, "deg")
        
          
        ### Save indices ###
        
        #Output directory
        outdir_mod <- paste0(outdir, "/baseline_", baseline_start[b], "_",  baseline_end[b],
                             "/", experiments[exp], "/", bc_methods[b], "/", models[m])
        
        dir.create(outdir_mod, recursive=TRUE)
        
        
        #Save as R object
        saveRDS(monthly_indices, paste0(outdir_mod, "/", outfile))
        
      }
      
    }
  }

}
