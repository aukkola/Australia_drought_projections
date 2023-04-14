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



#Set percentile and scale
percentile <- "Perc_15"

scale      <- 3

baseline   <- "1970-2005"

#Warming levels
envelopes =c(1,2)

#Variables
vars_all <- list(pr   = list(CMIP6 = "pr", BARPA = "pr", AWRA = "pr"),
                 qtot = list(CMIP6 = "mrro", BARPA = "mrro", AWRA = "qtot"),
                 sm   = list(CMIP6 = "mrsol_std_3m", BARPA = "mrso", AWRA = "sm"))
             
vars <- names(vars_all)

# var_labels <- c(pr   = "Precipitation", 
#                 qtot = "Runoff", 
#                 sm   = "Soil moisture") #labels for plotting


datasets <- c("CMIP6", "BARPA", "CABLE", "AWRA")


#List metrics
metrics <- c("duration", "rel_intensity", "time_in_drought")#, "frequency")


outdir <- paste0(path, "/Figures/")
dir.create(outdir)


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
                  time_in_drought = c(-100, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5,0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 100))


unit <- c(duration        = "months", 
          rel_intensity   = "% points",
          time_in_drought = "% points")



#Stippling settings
lwd <- 0.1
cex <- 0.15

panel_cex=0.7


#Need to resample CMIP6 models to the same resolution
#Create a 1 x 1 deg raster for this
res_raster <- raster(ext=extent(c(112.5, 153.75, -43.75, -10)), resolution=1)



#Loop through metrics
for (m in 1:length(metrics)) {
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_ensemble_changes_in_", metrics[m], "_",
               percentile, "_", scale, "_", vars[v], ".png"),
        height=7, width=8.3, units="in", res=400)
    
    
    par(mai=c(0.1, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.3, 0.4, 0.1))
    
    layout(matrix(c(1:3, 4, 4, 4, 5:7, 8, 8, 8), nrow=4, byrow=TRUE), heights=c(1, 0.3, 1, 0.3,
                                                                                1, 0.3, 1, 0.3,
                                                                                1, 0.3, 1, 0.3))
    
    
    #par(mfcol=c(3, 3))
    par(bty="n")
    
    
    #####################
    ### Read datasets ###
    #####################
 
    plot_data_hist <- list()
    plot_data_fut  <- list()

    
    for (d in 1:length(datasets)) {
      
      
      plot_data_hist[[datasets[d]]] <- brick()
      plot_data_fut[[datasets[d]]] <- lapply(envelopes, function(x) brick())

      
      #Find models (use historical to find them, only want to include models that have historical data)
      #should be all of them but just in case
      models <- list.files(paste0(path, "/Mean_drought_metrics_old_delete/", datasets[d], 
                                       "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                                       "/", metrics[m], "/historical/"))
      
      
      for (mod in 1:length(models)) {
        
        mod_path <-  paste0(path, "/Mean_drought_metrics_old_delete/", datasets[d], 
                            "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                            "/", metrics[m])
        


        
        
        #Have bc methods in addition to GCMs
        if (datasets[d] %in% c("AWRA", "CABLE")) {
          
          #Find BC methods
          bc_methods <- 
          
          
        ##################
        ### CMIP6 data ###
        ##################
          
        } else if (datasets[d] == "CMIP6") {
          
          #Find ensembles
          ensembles <- list.files(paste0(mod_path, "/historical/", models[mod]))
          
          #Initialise
          all_hist <- brick()
          all_fut  <- lapply(envelopes, function(x) brick())
          
          #Loop through ensembles
          for (ens in 1:length(ensembles)) {
            
            #Find historical file
            hist_file <- list.files(paste0(mod_path, "/historical/", models[mod], "/",
                                           ensembles[ens]), pattern=".nc", full.names=TRUE)
            
            #Should only find one file, check
            if (length(hist_file) != 1) stop("wrong CMIP6 historical file")
            
            #Read data
            all_hist <- addLayer(all_hist, raster(hist_file))
            
            
            ### Future data ##
            
            for(gw in 1:length(envelopes)) {
              
              #1 degree of warming
              files_deg <- list.files(paste0(mod_path, "/", envelopes[gw], "deg/", models[mod], "/", 
                                             ensembles[ens]),  recursive=TRUE, pattern=".nc", full.names=TRUE)
              
              #Check that found files, some GW levels don't have any
              if (length(files_deg) > 0) {
                
                #Read future data and calculate difference to historical
                data_deg <- brick(lapply(files_deg, raster)) - all_hist[[ens]]
                
                all_fut[[gw]] <- addLayer(all_fut[[gw]], calc(data_deg, median))
                
              }
               
            }
  
 
          } #ensembles
          
          #Calculate median of ensemble members and resample to a common resolution
          ens_median_hist <- resample(calc(all_hist, median), res_raster)
          ens_median_fut  <- lapply(all_fut, function(x) resample(calc(x, median), res_raster))
                                    
          
          plot_data_hist[[datasets[d]]] <- addLayer(plot_data_hist[[datasets[d]]], ens_median_hist)

          plot_data_fut[[datasets[d]]] <- lapply(envelopes, function(x) addLayer(plot_data_fut[[datasets[d]]][[x]],
                                                                                 ens_median_fut[[x]]))
          
        } #CMIP6
        
        
        
        ## BARPA ###
        

        
      }
      
      
      
      
      
      
   
      
      
      #Need to deal with multiple ensemble members
      
      
      
      
      
      
    }
    
    
    
    if (metrics[m] == "timing") {
      
      ############
      ### AWRA ###
      ############
      
      
      data_files_awra <- list.files(paste0(path, "/drought_metrics/", 
                                           scale, "-month/MRNBC/"),
                                    pattern=glob2rx(paste0("drought_metrics_*", vars[v], "*rcp45*.nc")),
                                    full.names=TRUE, recursive=TRUE)
      
      #CO2
      data_files_CO2  <- list.files(paste0(path, "/drought_metrics_CABLE/", 
                                           scale, "-month/CO2/"),
                                    pattern=glob2rx(paste0("drought_metrics_CABLE_CO2_*", vars[v], 
                                                           ".nc")),
                                    full.names=TRUE, recursive=TRUE)
      
      
      #noCO2
      data_files_noCO2  <- list.files(paste0(path, "/drought_metrics_CABLE/", 
                                             scale, "-month/noCO2/"),
                                      pattern=glob2rx(paste0("drought_metrics_CABLE_noCO2_*", vars[v], 
                                                             ".nc")),
                                      full.names=TRUE, recursive=TRUE)
      
      
      #Check number of files
      if (any (c(length(data_files_noCO2), length(data_files_CO2), length(data_files_awra)) != 4)){
        stop("Wrong number of files")
      }
      
      
      #Load data
      data_awra   <- lapply(data_files_awra, brick, varname=metrics[m])
      data_CO2    <- lapply(data_files_CO2, brick, varname=metrics[m])
      data_noCO2  <- lapply(data_files_noCO2, brick, varname=metrics[m])
      
      #Calculate historical and future mean
      
      #First get date indices so that use 1970-2005 for historical
      #and 2064-2099 for future
      
      dates <- getZ(data_awra[[1]])
      
      start_hist <- which(dates  == "1970-01-16")
      end_hist <- which(dates  == "2005-12-16")
      
      len_time <- length(start_hist:end_hist)
      
      start_rcp45 <- which(dates  == "2064-01-16")
      end_rcp45 <- which(dates  == "2099-12-16")
      
      
      #Historical mean (percentage of time in drought)
      data_hist_awra   <- brick(lapply(data_awra, function(x) sum(x[[start_hist:end_hist]], na.rm=TRUE) / len_time * 100))
      data_hist_CO2    <- brick(lapply(data_CO2, function(x) sum(x[[start_hist:end_hist]], na.rm=TRUE) / len_time * 100))
      data_hist_noCO2  <- brick(lapply(data_noCO2, function(x) sum(x[[start_hist:end_hist]], na.rm=TRUE) / len_time * 100))
      
      #Future mean
      data_rcp45_awra  <- brick(lapply(data_awra, function(x) sum(x[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time * 100))
      data_rcp45_CO2   <- brick(lapply(data_CO2, function(x) sum(x[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time * 100))
      data_rcp45_noCO2 <- brick(lapply(data_noCO2, function(x) sum(x[[start_rcp45:end_rcp45]], na.rm=TRUE) / len_time * 100))
      
      
    } else {
      
      
      ############
      ### AWRA ###
      ############
      
      
      data_files_hist_awra <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                                scale, "/historical/", vars[v], "/MRNBC/"),
                                         pattern=paste0("Mean_", metrics[m]),
                                         full.names=TRUE, recursive=TRUE)
      
      #Read datasets
      data_files_rcp45_awra <- list.files(paste0(path, "/Mean_drought_metrics/scale_", 
                                                 scale, "/rcp45/", vars[v],"/MRNBC/"),
                                          pattern=paste0("Mean_", metrics[m]),
                                          full.names=TRUE, recursive=TRUE)
      
      #Check number of files
      if (any (c(length(data_files_hist_awra), length(data_files_rcp45_awra)) != 4)){
        stop("Wrong number of CO2 files")
      }
      
      #CO2
      data_files_hist_CO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                scale, "/", vars[v], "/CO2/"),
                                         pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                "*_CO2_*_historical_*.nc")),
                                         full.names=TRUE, recursive=TRUE)
      
      data_files_rcp45_CO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                scale, "/", vars[v], "/CO2/"),
                                         pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                "*_CO2_*_rcp45_*.nc")),
                                         full.names=TRUE, recursive=TRUE)
      
      #Check number of files
      if (any (c(length(data_files_hist_CO2), length(data_files_rcp45_CO2)) != 4)){
        stop("Wrong number of CO2 files")
      }
      
      
      #noCO2
      data_files_hist_noCO2  <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                  scale, "/", vars[v], "/noCO2/"),
                                           pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                  "*_noCO2_*_historical_*.nc")),
                                           full.names=TRUE, recursive=TRUE)
      
      data_files_rcp45_noCO2 <- list.files(paste0(path, "/Mean_drought_metrics_CABLE/scale_", 
                                                  scale, "/", vars[v], "/noCO2/"),
                                           pattern=glob2rx(paste0("Mean_CABLE_", metrics[m], 
                                                                  "*_noCO2_*_rcp45_*.nc")),
                                           full.names=TRUE, recursive=TRUE)
      
      #Check number of files
      if (any (c(length(data_files_hist_noCO2), length(data_files_rcp45_noCO2)) != 4)){
        stop("Wrong number of CO2 files")
      }
      
      
      #Load data
      data_hist_awra    <- brick(lapply(data_files_hist_awra, raster))
      data_hist_CO2     <- brick(lapply(data_files_hist_CO2, raster))
      data_hist_noCO2   <- brick(lapply(data_files_hist_noCO2, raster))
      
      data_rcp45_awra   <- brick(lapply(data_files_rcp45_awra, raster))
      data_rcp45_CO2    <- brick(lapply(data_files_rcp45_CO2, raster))
      data_rcp45_noCO2  <- brick(lapply(data_files_rcp45_noCO2, raster))
      
    }
    
    
    
    #Calculate future change
    future_diff_rcp45_awra   <- data_rcp45_awra - data_hist_awra
    future_diff_rcp45_CO2    <- data_rcp45_CO2 - data_hist_CO2
    future_diff_rcp45_noCO2  <- data_rcp45_noCO2 - data_hist_noCO2
    
    #Calculate ensemble median change
    ens_median_awra  <- calc(future_diff_rcp45_awra, median)
    ens_median_CO2   <- calc(future_diff_rcp45_CO2, median)
    ens_median_noCO2 <- calc(future_diff_rcp45_noCO2, median)
    
    
    
    ################
    ### Plotting ###
    ################
    
    
    ### Ensemble median future change ###
    
    plot_data <- list(awra  = ens_median_awra,
                      CO2   = ens_median_CO2,
                      noCO2 = ens_median_noCO2)
    
    
    labs_mod <- c("AWRA", "CABLE CO2", "CABLE noCO2")   
    
    
    #Loop through experiments
    for (p in 1:length(plot_data)) {
      
      
      #Plot
      len <- length(lims_diff[[metrics[m]]])
      image(plot_data[[p]], breaks=lims_diff[[metrics[m]]], 
            col=c(cols_diff(len-1)),
            axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #Model label
      mtext(side=3, line=0, font=2, text=labs_mod[p], xpd=NA)
      
      
      if (p==1) mtext(side=2, line=1, text="Ensemble median")
      
    }
    
    
    
    #Empty plot
     
    plot(1, type="n", bty="n", yaxt="n", xaxt="n")
    
    #Legend
    len1 <- length(lims_diff[[metrics[m]]])-1
    add_raster_legend2(cols=cols_diff(len1), limits=lims_diff[[metrics[m]]][2:len1],
                       main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                       title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
  
    

    
    ### Model differences ###

    plot_data <- list(awra_CO2   = ens_median_awra - ens_median_CO2, #(ens_median_awra - ens_median_CO2) / ens_median_CO2 *100,
                      awra_noCO2 = ens_median_awra - ens_median_noCO2,
                      CO2_noCO2  = ens_median_noCO2 - ens_median_CO2)
    
    
    labs <- c("AWRA - CO2", "AWRA - noCO2", "noCO2 - CO2")
    
    cols <- colorRampPalette(rev(c("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", 
                               "#c7eae5", "#80cdc1", "#35978f", "#01665e")))
    
    #Loop through experiments
    for (p in 1:length(plot_data)) {
      
      
      #Plot
      len <- length(lims_diff[[metrics[m]]])
      image(plot_data[[p]], breaks=lims_diff[[metrics[m]]], 
            col=c(cols(len-1)),
            axes=FALSE, ann=FALSE, asp=1)
      
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #Model label
      mtext(side=3, line=0, font=2, text=labs[p], xpd=NA)
      
      
      if (p==1) mtext(side=2, line=1, text="Difference")
      
      
    }
    
    
    ### Legend ###
    
    #Empty plot
    plot(1, type="n", bty="n", yaxt="n", xaxt="n")
    
    #Legend
    len1 <- length(lims_diff[[metrics[m]]])-1
    add_raster_legend2(cols=cols(len1), limits=lims_diff[[metrics[m]]][2:len1],
                       main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
                       title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
  


    dev.off ()
    
  } #metrics
  
} #variables






