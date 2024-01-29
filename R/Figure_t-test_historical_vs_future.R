library(raster)
library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"


#Source functions
source(paste0(path,"/scripts/R/functions/add_raster_legend.R"))
source(paste0(path,"/scripts/R/functions/mask_data.R"))
source(paste0(path,"/scripts/R/functions/wtd_ttest.R"))


#Get world shapefile (for masking) and crop to Australia
data("wrld_simpl", package="maptools")
Aus <- crop(wrld_simpl, extent(112.5, 153.75, -43.75, -10))


#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

baseline   <- "1970-2005"

#Warming levels
envelopes =c(1,2)

#Variables
vars_all <- list(pr   = list(CMIP6 = "pr", BARPA = "pr", AWRA = "pr", CORDEX="pr"),
                 qtot = list(CMIP6 = "mrro", BARPA = "mrro", AWRA = "qtot", CORDEX="mrro"),
                 sm   = list(CMIP6 = "mrsol_std_2.89m", BARPA = "mrso", AWRA = "sm", CORDEX="mrso"))

vars <- names(vars_all)

var_labels <- c(pr   = "Precipitation",
                qtot = "Runoff",
                sm   = "Soil moisture") #labels for plotting


datasets <- c("CMIP6", "BARPA", "CORDEX", "AWRA" )


#List metrics
metrics <- c("time_in_drought", "duration", "rel_intensity")#, "frequency")


outdir <- paste0(path, "/Figures/")
dir.create(outdir)


#Parallelise
#beginCluster(12)


#####################
### Plot settings ###
#####################


#Set plot colours

#Historical mean
cols_hist <- colorRampPalette(c("#ffffe5", "#fee391",
                                         "#fe9929", "#cc4c02"))
                                         

#Future difference
col_opts <- rev(brewer.pal(11, 'RdYlBu'))
#col_opts[6] <- "grey80"
cols_diff <- colorRampPalette(col_opts) 

#Limits
lims_hist <- list(pr =   list(duration  = c(1, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              time_in_drought = c(5, 10, 12, 14, 16, 18, 1000)),
                  qtot = list(duration  = c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              time_in_drought = c(5, 10, 12, 14, 16, 18, 1000)),
                  sm =   list(duration  = c(2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 1000),
                              rel_intensity = c(0, 30, 40, 50, 60, 70, 80, 90, 100),
                              time_in_drought = c(5, 10, 12, 14, 16, 18, 1000)))

lims_diff <- list(duration        = c(-1000, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 1000),
                  rel_intensity   = c(-10000, -8, -6, -4, -2, 0,  2, 4, 6, 8, 10000),
                  time_in_drought = c(-100, -10, -8, -6 , -4, -2, 0 , 2, 4, 6, 8, 10, 100))


unit <- c(duration        = "months", 
          rel_intensity   = "% points",
          time_in_drought = "% points")



#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7


#Need to resample CMIP6 models to the same resolution
#Create a 1 x 1 deg raster for this
#res_raster <- raster(ext=extent(c(112.5, 153.75, -43.75, -10)), resolution=1)

#Then need to resample everything to AWRA/CABLE 5km resolution, create a raster for this
res_raster <- raster(ext=extent(c(112.5, 153.75, -43.75, -10)), resolution=0.05)


#Create folder to save plot data
plot_data_dir <- paste0(path, "/Plot_data/")


#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  # #One fewer variable for precip, set lower plot height    
  # if (vars[v] == "pr") {
  #   height=8
  # } else {
  height=6
  #}
  
  ### Set up figure ###
  png(paste0(outdir, "/FigureX", "_weighted_mean_all_datasets_in_", metrics[m], "_",
             percentile, "_", scale, "_all_vars_with_t-test.png"),
      height=height, width=8.3, units="in", res=400)
  
  
  par(mai=c(0.1, 0.1, 0.2, 0.1))
  par(omi=c(0.1, 0.3, 0.4, 0.1))
  
  
  #One fewer dataset for precip
  # if (vars[v] == "pr") {
  #   layout(matrix(c(1:9, 10, 11, 11), nrow=4, byrow=TRUE), heights=c(rep(1,3), 0.3) )
  # } else {
  #layout(matrix(c(1:3, 4,4,4, 5:7, 8,8,8), nrow=4, byrow=TRUE), heights=c(1, 0.3, 1, 0.3) )
  layout(matrix(c(1, 3, 5, 2, 4,6, 7,7,7), nrow=3, byrow=TRUE), heights=c(1, 1, 0.3) )
  
  #    }
  
  
  #par(mfrow=c(3, 3))
  par(bty="n")
  
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    #Progress
    print(paste0("Plotting variable: ", vars[v], ", metric ", m, "/", length(metrics)))
    
    
    #Save data for faster processing
    dataset_dir <- paste0(plot_data_dir, "/Combined/", vars[v], "/", metrics[m])
    dir.create(dataset_dir, recursive=TRUE)
    
    
    
    #####################
    ### Read datasets ###
    #####################
    
    plot_data_all_hist     <- list()
    plot_data_all_fut      <- list()
    plot_data_all_fut_diff <- list()
    
    
    #Loop through datasets    
    for (d in 1:length(datasets)) {
      
      print(paste0("plotting dataset ", d, "/", length(datasets)))
      
      #Save data for faster processing
      hist_out_file <- paste0(dataset_dir, "/historical_data_", datasets[d], "_",
                              vars[v], "_", metrics[m], ".nc")
      
      fut_out_files <- lapply(envelopes,  function(x) paste0(dataset_dir, "/", x, "deg_data_", datasets[d], "_",
                                                             vars[v], "_", metrics[m], ".nc"))
      
      fut_diff_out_files <- lapply(envelopes,  function(x) paste0(dataset_dir, "/", x, "deg_data_diff_", datasets[d], "_",
                                                                  vars[v], "_", metrics[m], ".nc"))
      
      
      if (any(!file.exists(c(hist_out_file, unlist(fut_out_files), unlist(fut_diff_out_files))))) {
      
        #Skip CABLE if rainfall
        if (vars[v] == "pr" & datasets[d] == "CABLE") next
        
        
        #Initialise
        plot_data_hist      <- brick()
        plot_data_fut       <- lapply(envelopes, function(x) brick())
        plot_data_fut_diff  <- lapply(envelopes, function(x) brick())
        
        
        #Find models (use historical to find them, only want to include models that have historical data)
        #should be all of them but just in case
        #For AWRA/CABLE this finds the bc methods, doesn't matter for processing
        model_path <- paste0(path, "/Mean_drought_metrics/", datasets[d], 
                             "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                             "/", metrics[m], "/historical/")
        
        pattern <- ""
        if (datasets[d] == "CABLE") pattern="/CO2/" 
        
        models <- list.files(paste0(model_path, pattern))
        
      
        #model name files
        mod_name_hist_file <-paste0(dataset_dir, "/historical_mod_names_", datasets[d], "_",
                                    vars[v], "_", metrics[m], ".csv")
        mod_name_fut_file <- lapply(envelopes,  function(x) paste0(dataset_dir, "/", x, "deg_mod_names_", datasets[d], "_",
                                                                   vars[v], "_", metrics[m], ".csv"))
        
        
        #Collate model names 
        model_names_hist <- vector()
        model_names_fut <- lapply(envelopes, function(x) vector())
        
        
    
        for (mod in 1:length(models)) {
          
          mod_path <-  paste0(path, "/Mean_drought_metrics/", datasets[d], 
                              "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                              "/", metrics[m])
          
          
          
          ##################
          ### AWRA/CABLE ###
          ##################
          
          
          #Have bc methods/RCMs in addition to GCMs
          if (datasets[d] %in% c("AWRA", "CABLE", "CORDEX")) {
            
            
            #models correspond to BC methods here. Find GCMs to loop through
            #For CABLE finds CO2/noCO2, set to CO2 for this plotting
            if (datasets[d] == "CABLE") {
              gcms <- "CO2"
            } else {
              gcms <- list.files(paste0(mod_path, "/historical/", models[mod]))
            }
            
            
            #Initialise
            all_fut      <- lapply(envelopes, function(x) brick()) #future mean
            all_fut_diff <- lapply(envelopes, function(x) brick()) #future difference
            
            #Loop through ensembles
            for (gc in 1:length(gcms)) {
              
              #Find historical file
              #Folders in a different order, deal with this
              if (datasets[d] == "CABLE") {
                hist_file <- list.files(paste0(mod_path, "/historical/", gcms[gc], "/",
                                               models[mod]), pattern=".nc", full.names=TRUE)
                
              } else if (datasets[d] == "CORDEX") {
                hist_file <- list.files(paste0(mod_path, "/historical/", models[mod], "/",
                                               gcms[gc], "/r1i1p1/"), pattern=".nc", full.names=TRUE)
                
              } else {
                hist_file <- list.files(paste0(mod_path, "/historical/", models[mod], "/",
                                               gcms[gc]), pattern=".nc", full.names=TRUE)
              }
              
              
              #Should only find one file, check
              if (length(hist_file) != 1) stop("wrong CABLE/AWRA historical file")
              
              #Read data
              hist_data <- raster(hist_file)
              
              
              ### Future data ##
              
              for(gw in 1:length(envelopes)) {
                
                #1 degree of warming
                
                
                #Folders in a different order, deal with this
                if (datasets[d] == "CABLE") {
                  files_deg <- list.files(paste0(mod_path, "/", envelopes[gw], "deg/", gcms[gc], "/", 
                                                 models[mod]),  recursive=TRUE, pattern=".nc", full.names=TRUE)
                  
                } else {
                  files_deg <- list.files(paste0(mod_path, "/", envelopes[gw], "deg/", models[mod], "/", 
                                                 gcms[gc]),  recursive=TRUE, pattern=".nc", full.names=TRUE)
                }
                
                
                
                
                #Check that found files, some GW levels don't have any
                if (length(files_deg) > 0) {
                  
                  #Read future data
                  data_deg <- brick(lapply(files_deg, raster)) 
                  
                  #calculate difference to historical
                  data_deg_diff <- data_deg - hist_data
                  
                  #calculate median
                  all_fut[[gw]] <- calc(data_deg, median) #future mean
                  
                  all_fut_diff[[gw]] <- calc(data_deg_diff, median) #future difference
                  
                  #save model names
                  model_names_fut[[gw]]  <- append(model_names_fut[[gw]], paste0(models[mod], "/", gcms[gc]))
                  
                }
                
              }
              
              
              #save model names
              model_names_hist <- append(model_names_hist, paste0(models[mod], "/", gcms[gc]))
              
              #Need to resample CORDEX to 5km resolution
              if (datasets[d] == "CORDEX") {
                
                #Calculate factor for disaggregating
                data_res <- res(hist_data) #only use first list element (GW envelope), they all match
                
                factor <- min(floor(data_res/ res(res_raster)))
                
                #First disaggregate coarser data to a resolution as close to 5km as possible.
                #Then resample using bilinear. This will hopefully minimise the creation of new data
                #by repeating coarse data for finer pixels
                hist_data    <- resample(disaggregate(hist_data, fact=factor), res_raster)
                all_fut      <- lapply(all_fut, function(x) resample(disaggregate(x, fact=factor), res_raster))
                all_fut_diff <- lapply(all_fut_diff, function(x) resample(disaggregate(x, fact=factor), res_raster))
                
              }
              
              
              
              #Collate
              #historical mean
              plot_data_hist <- addLayer(plot_data_hist, hist_data)
          
              #future mean
              plot_data_fut <- lapply(envelopes, function(x) if(!is.null(all_fut[[x]]))
                                                                    addLayer(plot_data_fut[[x]],
                                                                             all_fut[[x]])
                                                                    else plot_data_fut[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
                      
              #future difference
              plot_data_fut_diff <- lapply(envelopes, function(x) if(!is.null(all_fut_diff[[x]]))
                                                                  addLayer(plot_data_fut_diff[[x]],
                                                                           all_fut_diff[[x]])
                                                                  else plot_data_fut_diff[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
                                                                
            } #GCMs            
            
            
            
            
            ##################
            ### CMIP6 data ###
            ##################
            
          } else if (datasets[d] %in% c("CMIP6", "BARPA")) {
            
            #Find ensembles
            ensembles <- list.files(paste0(mod_path, "/historical/", models[mod]))
            
            #Initialise
            all_hist      <- brick()
            all_fut       <- lapply(envelopes, function(x) brick())
            all_fut_diff  <- lapply(envelopes, function(x) brick())
            
            #Loop through ensembles
            for (ens in 1:length(ensembles)) {
              
              #Find historical file
              hist_file <- list.files(paste0(mod_path, "/historical/", models[mod], "/",
                                             ensembles[ens]), pattern=".nc", full.names=TRUE)
              
              #Should only find one file, check
              if (length(hist_file) < 1) stop("wrong CMIP6 historical file")
              
              #Bug in processing code, accidentally wrote multiple historical files
              #take the first one if found multiple
              if (length(hist_file) > 1) hist_file <- hist_file[1]
              
              #Read data
              all_hist <- addLayer(all_hist, raster(hist_file))
              
              
              ### Future data ##
              
              for(gw in 1:length(envelopes)) {
                
                #1 degree of warming
                files_deg <- list.files(paste0(mod_path, "/", envelopes[gw], "deg/", models[mod], "/", 
                                               ensembles[ens]),  recursive=TRUE, pattern=".nc", full.names=TRUE)
                
                #Check that found files, some GW levels don't have any
                if (length(files_deg) > 0) {
  
                  #Read future data
                  data_deg <- brick(lapply(files_deg, raster)) 
                  
                  #calculate difference to historical
                  data_deg_diff <- data_deg - all_hist[[ens]]
                  
                  
                  #calculate median
                  all_fut[[gw]] <- calc(data_deg, median) #future mean
                  
                  all_fut_diff[[gw]] <- addLayer(all_fut_diff[[gw]], calc(data_deg_diff, median)) #!!!!!!!!!!!!!!! fix all of this
                  
                }
                
              }
              
              
            } #ensembles
            
            #Calculate median of ensemble members and resample to a common resolution
            ens_median_hist      <- resample(calc(all_hist, median), res_raster)
            ens_median_fut       <- lapply(all_fut, function(x) if(!all(is.na(values(x)))) resample(calc(x, median), res_raster))
            ens_median_fut_diff  <- lapply(all_fut_diff, function(x) if(!all(is.na(values(x)))) resample(calc(x, median), res_raster))
            
            
            ### Resample ###
            
            #Calculate factor for disaggregating
            data_res <- res(ens_median_hist) #only use first list element (GW envelope), they all match
            
            factor <- min(floor(data_res/ res(res_raster)))
            
            #First disaggregate coarser data to a resolution as close to 5km as possible.
            #Then resample using bilinear. This will hopefully minimise the creation of new data
            #by repeating coarse data for finer pixels
            ens_median_hist     <- resample(disaggregate(ens_median_hist, fact=factor), res_raster)
            ens_median_fut      <- lapply(ens_median_fut, function(x) if(!is.null(x))
                                                                           resample(disaggregate(x, fact=factor), res_raster)
                                                                           else x)
            ens_median_fut_diff <- lapply(ens_median_fut_diff, function(x) if(!is.null(x))
                                                                           resample(disaggregate(x, fact=factor), res_raster)
                                                                           else x)
            
            
            ### Collate ###
            
            #historical mean
            plot_data_hist <- addLayer(plot_data_hist, ens_median_hist)
            
            #future mean
            plot_data_fut <- lapply(envelopes, function(x) if(!is.null(ens_median_fut[[x]]))
                                                            addLayer(plot_data_fut[[x]],
                                                                     ens_median_fut[[x]])
                                                            else plot_data_fut[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
            #future difference                                              
            plot_data_fut_diff <- lapply(envelopes, function(x) if(!is.null(ens_median_fut_diff[[x]]))
                                                                  addLayer(plot_data_fut_diff[[x]],
                                                                           ens_median_fut_diff[[x]])
                                                                  else plot_data_fut_diff[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
                                                                
            #save model names
            model_names_hist <- append(model_names_hist, models[mod])
            
            for (ii in 1:length(envelopes)) {
              if(!is.null(ens_median_fut[[ii]])) {
                model_names_fut[[ii]]  <- append(model_names_fut[[ii]], models[mod])
              }
            }
            
            
          } #CMIP6/BARPA
        } #models
        
        
        #Add layer names
        names(plot_data_hist) <- model_names_hist
        
        for (mm in 1:length(plot_data_fut_diff))  {
          names(plot_data_fut[[mm]]) <- model_names_fut[[mm]]
          names(plot_data_fut_diff[[mm]]) <- model_names_fut[[mm]]
  
        }
  
  # 
  # 
  #         #Write outputs
  #         writeRaster(plot_data_hist, hist_out_file, overwrite=TRUE, varname=metrics[m])
  #         lapply(1:length(plot_data_fut_diff), function(x) writeRaster(plot_data_fut_diff[[x]], fut_out_files[[x]],
  #                                                                      overwrite=TRUE, varname=metrics[m]))
  # 
  # 
  #         #Also need to write out model names as these don't get saved into the file otherwise
  #         write.csv(model_names_hist, mod_name_hist_file)
  # 
  #         lapply(1:length(model_names_fut), function(x) write.csv(model_names_fut[[x]], mod_name_fut_file[[x]]))
  # 
  # 
  #         #If files already processed, read them in
  #       } else {
  # 
  #         #Read data in
  #         plot_data_hist <- brick(hist_out_file)
  #         plot_data_fut_diff  <- lapply(fut_out_files, brick)
  # 
  #         #Also get model names (need to grab second column as first one is just indices)
  #         model_names_hist  <- read.csv(mod_name_hist_file)[,2]
  #         model_names_fut   <- lapply(mod_name_fut_file, function(x) read.csv(x)[,2])
  # 
  #       }
  # 
      
        
        ### Find common models
        
        #Match models for historical and future runs
        #First collate historical and future names
        all_names <- model_names_fut
        all_names[[length(all_names)+1]] <- model_names_hist
        
        #Then find common models
        common_models <- Reduce(intersect, all_names) #needs capital R, found on stack overflow
        
    
        #Then collate datasets (selecting only common models)
        plot_data_all_hist[[d]]     <- plot_data_hist[[which(model_names_hist %in% common_models)]]
        plot_data_all_fut[[d]]      <- lapply(envelopes, function(env) plot_data_fut[[env]][[which(model_names_fut[[env]] %in% common_models)]])
        plot_data_all_fut_diff[[d]] <- lapply(envelopes, function(env) plot_data_fut_diff[[env]][[which(model_names_fut[[env]] %in% common_models)]])
        
        
        
        #Write outputs
        
        #historical
        writeRaster(plot_data_all_hist[[d]], hist_out_file, overwrite=TRUE, varname=metrics[m])
        
        #future mean
        
        lapply(1:length(plot_data_all_fut[[d]]), function(x) writeRaster(plot_data_all_fut[[d]][[x]], fut_out_files[[x]], 
                                                                overwrite=TRUE, varname=metrics[m]))
        
        #future difference
        lapply(1:length(plot_data_all_fut_diff[[d]]), function(x) writeRaster(plot_data_all_fut_diff[[d]][[x]], fut_diff_out_files[[x]], 
                                                                         overwrite=TRUE, varname=metrics[m]))
          
 
      } else {
        
        #Read data in
        plot_data_all_hist[[d]] <- brick(hist_out_file)
        plot_data_all_fut[[d]] <- lapply(fut_out_files, brick)
        plot_data_all_fut_diff[[d]] <- lapply(fut_diff_out_files, brick)
        
        
      } #if file exists
    } #datasets
    
    
    
    ################
    ### Plotting ###
    ################
    
    #Calculate weights for all datasets
    
    #Plot global warming envelopes
    for(env in envelopes) {
      
      ### Calculate ensemble mean ###
      
      #Calculate weights (weights don't add up to 1 but this doesn't matter for the R function,
      #checked that returns the same results as normalising weights to equal 1)
      weights <- unlist(lapply(plot_data_all_fut_diff, function(x) rep(1/nlayers(x[[env]]), nlayers(x[[env]]))))
      
      #Get data for GW envelope (resample and crop to make extent same, thsi is because AWRA coordinates are slightly shifted)
      all_data_diff <- brick(lapply(plot_data_all_fut_diff, function(x) resample(crop(x[[env]], res_raster), res_raster)))
      
      #Calculate weighted ensemble mean
      ensemble_mean <- weighted.mean(all_data_diff, w=weights, na.rm=TRUE)
      
      
      ### Calculate t-test ###
      
      #Calculate t-test (test if historical and future means are significantly different)
      all_data_hist <- brick(lapply(plot_data_all_hist, function(x) resample(crop(x, res_raster), res_raster)))
      
      all_data_fut  <- brick(lapply(plot_data_all_fut, function(x) resample(crop(x[[env]], res_raster), res_raster)))
      
      combined_hist_fut <- brick(list(all_data_hist, all_data_fut))
      
      #First check that same amount of layers in hist and future data
      #(t-test function splits them and assumes the same number of layers)
      if(nlayers(all_data_hist) != nlayers(all_data_fut)) stop("different nlayers in hist and fut")
      
      #Calculate future test
      ttest_pvalue <- calc(combined_hist_fut, function(x) wtd_ttest_raster(x, weight=weights))
      
      
      # 
      # #Calculate model agreement (AR4 metric)
      # fut_mod_agr <- calc(all_data, function(x) perc_agree_on_sign_weighted_AR4_method(x, weights=weights))
      
      #Set pixels with no agreement to grey
      #First set them to a value smaller than min value in lims_diff
      mask_value <- min(lims_diff[[metrics[m]]]) -1
      ensemble_mean[ttest_pvalue > 0.05] <- mask_value
      
      #Set colours (grey for no agreement)
      len <- length(lims_diff[[metrics[m]]])
      
      plot_colours <- c("grey80", c(cols_diff(len-1)))
      
      
      #Plot
      image(mask(ensemble_mean, Aus), breaks=c(mask_value, lims_diff[[metrics[m]]]), 
            col=plot_colours,
            axes=FALSE, ann=FALSE, asp=1)
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #time period label
      if (v == 1) mtext(side=2, paste0(env, " degree"))
      
      if (env ==1) mtext(side=3, var_labels[vars[v]])
    }
    
    
    
    
    # #Again mask where more than 20% of pixels are missing
    # mask_data_fut <- calc(plot_data_fut_diff[[env]][[which(model_names_fut[[env]] %in% common_models)]], mask_data)
    # 
    # #calculate median
    # fut_median <- calc(mask_data_fut, fun=median, na.rm=TRUE)
    # 
    
    
    #Future legend
    
    if (v == length(vars)) {
      
      #Empty plot
      plot(1, type="n", bty="n", yaxt="n", xaxt="n")
      
      #Legend
      len1 <- length(lims_diff[[metrics[m]]])-1
      add_raster_legend2(cols=cols_diff(len1), limits=lims_diff[[metrics[m]]][2:len1],
                         main_title=unit[metrics[m]], plot_loc=c(0.2,0.8,0.63, 0.77),
                         title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
      
      
    }
    
    
  } #variables
  
  dev.off ()
  
} #metrics






