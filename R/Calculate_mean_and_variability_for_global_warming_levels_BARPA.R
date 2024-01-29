library(raster)
library(lubridate)

#clear R environment
rm(list=ls(all=TRUE))


path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"

#Source function
source(paste0(path, "/scripts/R/functions/freq.R"))


#Output directory
outdir <- paste0(path, "/Mean_and_variability/")

scale <- 3

warming_levels <- c(1,2)

#Historical period
hist_ref <- c(1970, 2005)

metrics <- list(mean="mean", 
                variability="variability", 
                CV="CV")

units <- list(mean="mm/month", 
              variability="mm/month", 
              CV="%")


#Parallelise
beginCluster(12)



#############
### BARPA ###
#############

print("Processing BARPA")

barpa_path <- paste0(path, "/Drought_metrics/BARPA/")

experiments <- list.files(barpa_path, pattern="ssp")


#Loop through experimens
for (exp in 1:length(experiments)) {
  
  variables <- list.files(paste0(barpa_path, "/", experiments[exp]))
  
  #Loop through variables
  for (v in 1:length(variables)) {
    
    dir <- paste0(barpa_path, "/", experiments[exp], "/", variables[v],
                  "/Perc_15/Baseline_1970_2005/Scale_", scale)
    
    gcms <- list.files(dir)
    
    #Loop through GCMs
    for (g in 1:length(gcms)) {
      
      ensembles <- list.files(paste0(dir, "/", gcms[g]))
      
      #loop through ensembles
      for (ens in 1:length(ensembles)) {
        
        #Find historical file
        hist_file <- list.files(paste0(barpa_path, "/historical/", variables[v],
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
          
          
          #output directory for metric
          outdir_mod <- paste0(outdir, "/BARPA/scale_", scale, "/", variables[v],
                               "/", names(metrics)[m], "/") 
          
          
          ##################
          ### Historical ###
          ##################
          
          #Load historical data
          hist_data <- brick(hist_file, varname=metrics[[m]])
          
          years <- year(getZ(hist_data))
          
          #Check that years are reasonable (GFDL starts 1861 so use that, others
          #should be 1850)
          if (min(years) != 1960) stop(paste0("wrong historical start year, file: ",
                                              hist_file))  
          
          #Historical output file
          outdir_hist <- paste0(outdir_mod, "/historical/", gcms[g], "/",
                                ensembles[ens])
          dir.create(outdir_hist, recursive=TRUE)
          
          hist_outfile <- paste0(outdir_hist, "/Historical_mean_", variables[v], "_",
                                 names(metrics)[m], "_", gcms[g], "_", ensembles[ens], ".nc")
          
          if (!file.exists(hist_outfile)) {
            
            
            #Extract historical mean period
            start_ind <- which(years == hist_ref[1])[1]
            end_ind   <- tail(which(years == hist_ref[2]), n=1)
            
            if (metrics[[m]] == "timing") {
              
              #Calculate time under drought (in %)
              hist_mean <- clusterR(hist_data[[start_ind:end_ind]], calc, args=list(fun=freq))
              
            } else {
              #Take mean of the metric
              hist_mean <- clusterR(hist_data[[start_ind : end_ind]], mean, args=list(na.rm=TRUE))
            }
            
            
            #Write output        
            writeRaster(hist_mean, hist_outfile, overwrite=TRUE, varname=names(metrics)[m],
                        longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])

          }
            
          
          ###################
          ### Future data ###
          ###################
          
          #Load historical data
          fut_data <- brick(fut_file, varname=metrics[[m]])
          
          years_fut <- year(getZ(fut_data))
          
          if (years_fut[1] != 2015) stop(paste0("wrong future start year, file: ",
                                                fut_file))  
          
          #Combine historical and future data (global warming level indices are
          #relative to the start of the historical period)
          all_data <- addLayer(hist_data, fut_data)
          
          
          
          #Read index file for global warming levels
          gw_indices <- readRDS(list.files(paste0(path, "/Global_warming_levels/",
                                                  "BARPA/baseline_1970_2005/", 
                                                  experiments[exp], "/", gcms[g], "/", 
                                                  ensembles[ens]), pattern=".rds", full.names=TRUE))
          
          #Loop through global warming levels
          for (w in 1:length(warming_levels)) {
            
            #Get indices for this warming level
            inds <- gw_indices[[paste0(w, "deg")]]
            
            #if no indices, skip
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
              
          } #GW levels
        } #metrics
      } #ensembles
    } #GCMs
  } #variables
} #experiments

endCluster()



