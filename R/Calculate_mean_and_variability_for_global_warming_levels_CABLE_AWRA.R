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

metrics <- list(duration="duration", 
                rel_intensity="rel_intensity", 
                time_in_drought="timing")

units <- list(duration="months", 
              rel_intensity="%", 
              time_in_drought="%")

#Parallelise
beginCluster(12)


############
### AWRA ###
############

awra_path <- paste0(path, "/../Steven_CABLE_runs/drought_metrics/")

experiments <- c("rcp45", "rcp85")

variables <- c("qtot", "pr", "sm")


#Loop through experimens
for (exp in 1:length(experiments)) {
  
  
  bc_methods <- list.files(paste0(awra_path, "/", scale, "-month/"))
  
  
  #loop through BC methods
  for (bc in 1:length(bc_methods)) {
  
    gcms <- list.files(paste0(awra_path, "/", scale, "-month/", bc_methods[bc]))
    
    #Loop through GCMs
    for (g in 1:length(gcms)) {
    
    
      #Loop through variables
      for (v in 1:length(variables)) {
        
        
        #Find file (combined historical and future)
        file <- list.files(paste0(awra_path, "/", scale, "-month/", bc_methods[bc],
                                  "/", gcms[g]), pattern=paste0(variables[v], "_", experiments[exp]),
                           full.names=TRUE)
                           
                          
        #Check that found one file each
        if(length(file) != 1) stop("wrong number of cmip5 files")
        
        
        #Loop through metrics
        for (m in 1:length(metrics)) {
          
          
          #output directory for metric
          outdir_mod <- paste0(outdir, "/AWRA/scale_", scale, "/", variables[v],
                               "/", names(metrics)[m], "/") 
          
          #Get data
          data <- brick(file, varname=metrics[[m]])
          
          
          
          ##################
          ### Historical ###
          ##################
          
          years <- year(getZ(data))
          
          #Check that years are reasonable (GFDL starts 1861 so use that, others
          #should be 1850)
          if (min(years) != 1960) stop(paste0("wrong historical start year, file: ",
                                              file))  
          
          
          #Extract historical mean period
          start_ind <- which(years == hist_ref[1])[1]
          end_ind   <- tail(which(years == hist_ref[2]), n=1)
          
          
          #Save historical data
          outdir_hist <- paste0(outdir_mod, "/historical/", bc_methods[bc], "/", gcms[g])
          
          dir.create(outdir_hist, recursive=TRUE)
          
          hist_outfile <- paste0(outdir_hist, "/Historical_mean_", variables[v], "_",
                                 names(metrics)[m], "_", bc_methods[bc], "_", gcms[g], ".nc")
          
          
          if (!file.exists(hist_outfile)) {
            
            if (metrics[[m]] == "timing") {
              
              #Calculate time under drought (in %)
              hist_mean <- clusterR(data[[start_ind:end_ind]], calc, args=list(fun=freq))
              
            } else {
              #Take mean of the metric
              hist_mean <- clusterR(data[[start_ind : end_ind]], mean, args=list(na.rm=TRUE)) #mean(data[[start_ind : end_ind]], na.rm=TRUE)
            }
            
            writeRaster(hist_mean, hist_outfile, overwrite=TRUE, varname=names(metrics)[m],
                        longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
            
            
          } else {
            hist_mean <- raster(hist_outfile)
          }
           
          ###################
          ### Future data ###
          ###################
          
          
          #Read index file for global warming levels
          gw_indices <- readRDS(list.files(paste0(path, "/Global_warming_levels/",
                                                  "CABLE_AWRA/baseline_1970_2005/", 
                                                  experiments[exp], "/", bc_methods[bc], "/",
                                                  gcms[g]), pattern=".rds", full.names=TRUE))
          
          #Loop through global warming levels
          for (w in 1:length(warming_levels)) {
            
            #Get indices for this warming level
            inds <- gw_indices[[paste0(w, "deg")]]
            
            if (length(inds) == 0) next
            
            
            if (metrics[[m]] == "timing") {
              
              #Calculate time under drought (in %)
              fut_mean <- clusterR(data[[inds]], calc, args=list(fun=freq))
              
            } else {
              #Take mean of the metric
              fut_mean <- clusterR(data[[inds]], mean, args=list(na.rm=TRUE))
            }
            
            
            #Save future data
            outdir_fut <- paste0(outdir_mod, "/", w, "deg/", bc_methods[bc], "/", gcms[g])
            dir.create(outdir_fut, recursive=TRUE)
            
            fut_outfile <- paste0(outdir_fut, "/", w, "deg_mean_", variables[v], "_",
                                  names(metrics)[m], "_", bc_methods[bc], "_", gcms[g], "_",
                                  experiments[exp], ".nc")
            
            writeRaster(fut_mean, fut_outfile, overwrite=TRUE, varname=names(metrics)[m],
                        longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
            
            
          } #GW levels
        } #metrics
      } #variables
    } #GCMs
  } #BC methods
} #experiments
        



#############
### CABLE ###
#############


cable_path <- paste0(path, "/../Steven_CABLE_runs/drought_metrics_CABLE/")

variables <- c("qtot", "sm")


#Loop through experimens
for (exp in 1:length(experiments)) {
  
  co2 <- list.files(paste0(cable_path, "/", scale, "-month/"))
  
  #loop through CO2/noCO2
  for (c in 1:length(co2)) {
    
    gcms <- list.files(paste0(cable_path, "/", scale, "-month/", co2[c]))
    
    #Skip this combo as no simulations available
    if (co2[c] == "noCO2" & experiments[exp] == "rcp85") next
    
    
    #Loop through GCMs
    for (g in 1:length(gcms)) {
      
      
      #Loop through variables
      for (v in 1:length(variables)) {
        
        
        #Find file (combined historical and future)
        file <- list.files(paste0(cable_path, "/", scale, "-month/", co2[c],
                                  "/", gcms[g]), pattern=paste0(variables[v], "_", experiments[exp]),
                           full.names=TRUE)
        
        
        #Check that found one file each
        if(length(file) != 1) stop("wrong number of CABLE files")
        
        #Progress
        print(paste0("file: ", file)) 
        
        
        #Loop through metrics
        for (m in 1:length(metrics)) {
          
          
          #output directory for metric
          outdir_mod <- paste0(outdir, "/CABLE/scale_", scale, "/", variables[v],
                               "/", names(metrics)[m], "/") 
          
          #Get data
          data <- brick(file, varname=metrics[[m]])
          
          
          
          ##################
          ### Historical ###
          ##################
          
          years <- year(getZ(data))
          
          #Check that years are reasonable (GFDL starts 1861 so use that, others
          #should be 1850)
          if (min(years) != 1960) stop(paste0("wrong historical start year, file: ",
                                              file))  
          
          
          #Extract historical mean period
          start_ind <- which(years == hist_ref[1])[1]
          end_ind   <- tail(which(years == hist_ref[2]), n=1)
          
          #Save historical data
          outdir_hist <- paste0(outdir_mod, "/historical/", co2[c], "/", gcms[g])
          
          dir.create(outdir_hist, recursive=TRUE)
          
          hist_outfile <- paste0(outdir_hist, "/Historical_mean_", variables[v], "_",
                                 names(metrics)[m], "_", co2[c], "_", gcms[g], ".nc")
          
          
          if (!file.exists(hist_outfile)) {
            
            if (metrics[[m]] == "timing") {
              
              #Calculate time under drought (in %)
              hist_mean <- clusterR(data[[start_ind:end_ind]], calc, args=list(fun=freq))
              
            } else {
              #Take mean of the metric
              hist_mean <- clusterR(data[[start_ind : end_ind]], mean, args=list(na.rm=TRUE))
            }
            
            
            
            writeRaster(hist_mean, hist_outfile, overwrite=TRUE, varname=names(metrics)[m],
                        longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
          } else {
            hist_mean <- raster(hist_outfile)
          }

          
          ###################
          ### Future data ###
          ###################
          
          
          #Read index file for global warming levels (use MRNBC bc method as that was CABLE input)
          gw_indices <- readRDS(list.files(paste0(path, "/Global_warming_levels/",
                                                  "CABLE_AWRA/baseline_1970_2005/", 
                                                  experiments[exp], "/MRNBC/",
                                                  gcms[g]), pattern=".rds", full.names=TRUE))
          
          #Loop through global warming levels
          for (w in 1:length(warming_levels)) {
            
            #Get indices for this warming level
            inds <- gw_indices[[paste0(w, "deg")]]
            
            if (length(inds) == 0) next
            
            
            #Save future data
            outdir_fut <- paste0(outdir_mod, "/", w, "deg/", co2[c], "/", gcms[g])
            dir.create(outdir_fut, recursive=TRUE)
            
            fut_outfile <- paste0(outdir_fut, "/", w, "deg_mean_", variables[v], "_",
                                  names(metrics)[m], "_", co2[c], "_", gcms[g], "_",
                                  experiments[exp], ".nc")
            
            if (!file.exists(fut_outfile)) {
              
              if (metrics[[m]] == "timing") {
                
                #Calculate time under drought (in %)
                fut_mean <- clusterR(data[[inds]], calc, args=list(fun=freq))
                
              } else {
                #Take mean of the metric
                fut_mean <- clusterR(data[[inds]], mean, args=list(na.rm=TRUE))
              }
              
              
              #Write output
              writeRaster(fut_mean, fut_outfile, overwrite=TRUE, varname=names(metrics)[m],
                          longname=paste0("mean ", names(metrics)[m]), varunit=units[[m]])
 
            }
            
               
          } #GW levels
        } #metrics
      } #variables
    } #GCMs
  } #BC methods
} #experiments


endCluster()


