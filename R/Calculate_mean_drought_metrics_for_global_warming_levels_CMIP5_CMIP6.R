library(raster)
library(lubridate)

#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"

#Source function
source(paste0(path, "/scripts/R/functions/freq.R"))

#Output directory
outdir <- paste0(path, "/Mean_drought_metrics/")

scale <- 3

warming_levels <- c(1,2)

#Historical period
hist_ref <- c(1970, 2005)

metrics <- list(duration="duration", 
                rel_intensity="rel_intensity", 
                time_in_drought="timing")

units <- list(duration="months", 
              rel_intensity="%", 
              time_in_drought="%")

#Parallelise
beginCluster(12)


#######################
### CMIP5 and CMIP6 ###
#######################

dataset <- c("CMIP5", "CMIP6")
pattern <- c("rcp", "ssp")

#For checking that data starts in the right year
fut_start_yr=list(CMIP5=2006,
                  CMIP6=2015)

for (d in 1:length(dataset)) {
  
  #Progress
  print(paste0("Processing ", dataset[d]))
  
  cmip_path <- paste0(path, "/Drought_metrics/", dataset[d])
  
  experiments <- list.files(cmip_path, pattern=pattern[d])
  
  
  #Loop through experimens
  for (exp in 1:length(experiments)) {
    
    variables <- list.files(paste0(cmip_path, "/", experiments[exp]))
    
    #Loop through variables
    for (v in 1:length(variables)) {
      
      dir <- paste0(cmip_path, "/", experiments[exp], "/", variables[v],
                    "/Perc_15/Baseline_1970_2005/Scale_", scale)
      
      gcms <- list.files(dir)
      
      #Loop through GCMs
      for (g in 1:length(gcms)) {
        
        ensembles <- list.files(paste0(dir, "/", gcms[g]))
        
        #loop through ensembles
        for (ens in 1:length(ensembles)) {
          
          #Find historical file
          hist_file <- list.files(paste0(cmip_path, "/historical/", variables[v],
                                         "/Perc_15/Baseline_1970_2005/Scale_", scale, "/",
                                         gcms[g], "/", ensembles[ens]),
                                  pattern=".nc", full.names=TRUE)
          
          
          #Find future file
          fut_file <- list.files(paste0(dir, "/", gcms[g], "/", ensembles[ens]),
                                 pattern=".nc", full.names=TRUE)
          
          #Check that found one file each
          if(length(fut_file) != 1 | length(hist_file) != 1) stop("wrong number of cmip5 files")
          
          
          #Loop through metrics
          for (m in 1:length(metrics)) {
            
            start=Sys.time()
            
            #output directory for metric
            outdir_mod <- paste0(outdir, "/", dataset[d], "/scale_", scale, "/", variables[v],
                                 "/", names(metrics)[m], "/") 
            
            
            
            #First of all check that have global warming levels available. Some
            #ensemble members don't seem to have tas data available
            gw_file <- list.files(paste0(path, "/Global_warming_levels/",
                                         dataset[d], "/baseline_1970_2005/", 
                                         experiments[exp], "/", gcms[g], "/", 
                                         ensembles[ens]), pattern=".rds", full.names=TRUE)
            
            if(length(gw_file) == 0) {
              print(paste0("Skipping ", outdir_mod, ", no GW levels available"))
              next
            }
            
            #Read index file for global warming levels
            gw_indices <- readRDS(gw_file)
            
            
            ##################
            ### Historical ###
            ##################
            
            #Load historical data
            hist_data <- brick(hist_file, varname=metrics[[m]])
            
            years <- year(getZ(hist_data))
            
            #Check that years are reasonable (GFDL starts 1861 so use that, others
            #should be 1850)
            if (min(years) > 1861) stop(paste0("wrong historical start year, file: ",
                                               hist_file))  
            
            
            #Historical output file
            outdir_hist <- paste0(outdir_mod, "/historical/", gcms[g], "/",
                                  ensembles[ens])
            dir.create(outdir_hist, recursive=TRUE)
            
            hist_outfile <- paste0(outdir_hist, "/Historical_mean_", variables[v], "_",
                                   names(metrics)[m], "_", gcms[g], "_", ensembles[ens], ".nc")
            
            print(paste0("processing: ", hist_outfile))
            
            if (!file.exists(hist_outfile)) {
              
              #Extract historical mean period
              start_ind <- which(years == hist_ref[1])[1]
              end_ind   <- tail(which(years == hist_ref[2]), n=1)
              
              if (metrics[[m]] == "timing") {
                
                #Calculate time under drought (in %)
                hist_mean <- clusterR(hist_data[[start_ind:end_ind]], calc, args=list(fun=freq))
                
              } else {
                #Take mean of the metric
                hist_mean <- clusterR(hist_data[[start_ind : end_ind]], mean, args=list(na.rm=TRUE)) #mean(data[[start_ind : end_ind]], na.rm=TRUE)
              }
              
              
              #Write output            
              writeRaster(hist_mean, hist_outfile, overwrite=TRUE, varname=names(metrics)[m],
                          longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
              
            } else {
              hist_mean <- raster(hist_outfile)
            }
              
            
            ###################
            ### Future data ###
            ###################
            
            #Load historical data
            fut_data <- brick(fut_file, varname=metrics[[m]])
            
            years_fut <- year(getZ(fut_data))
            
            
            if (years_fut[1] != fut_start_yr[dataset[d]]) stop(paste0("wrong future start year, file: ",
                                                  fut_file))  
            
            #Combine historical and future data (global warming level indices are
            #relative to the start of the historical period)
            all_data <- addLayer(hist_data, fut_data)
            
            
            
             #Loop through global warming levels
            for (w in 1:length(warming_levels)) {
              
              #Get indices for this warming level
              inds <- gw_indices[[paste0(w, "deg")]]
              
              if (length(inds) == 0) next
              
              #Future output file
              outdir_fut <- paste0(outdir_mod, "/", w, "deg/", gcms[g], "/",
                                   ensembles[ens])
              dir.create(outdir_fut, recursive=TRUE)
              
              fut_outfile <- paste0(outdir_fut, "/", w, "deg_mean_", variables[v], "_",
                                    names(metrics)[m], "_", gcms[g], "_", ensembles[ens], "_",
                                    experiments[exp], ".nc")
              
              if (!file.exists(fut_outfile)) {
                
                if (metrics[[m]] == "timing") {
                  
                  #Calculate time under drought (in %)
                  fut_mean <- clusterR(all_data[[inds]], calc, args=list(fun=freq))
                  
                } else {
                  #Take mean of the metric
                  fut_mean <- clusterR(all_data[[inds]], mean, args=list(na.rm=TRUE))
                }
                
                
                #Write output
                writeRaster(fut_mean, fut_outfile, overwrite=TRUE, varname=names(metrics)[m],
                            longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
                
                
              }
              
              end=Sys.time()
              print(end-start)
              
            } #GW levels
          } #metrics
        } #ensembles
      } #GCMs
    } #variables
  } #experiments
} #datasets

endCluster()
