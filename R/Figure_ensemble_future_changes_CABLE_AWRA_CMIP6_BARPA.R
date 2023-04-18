library(raster)
library(RColorBrewer)
library(maptools)
library(maps)

#clear R environment
rm(list=ls(all=TRUE))

## NEED TO CHANGE FILE PATH
#AND ADD CHECK THAT SAME MODELS INCLUDED IN HIST AND FUTURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


datasets <- c("CMIP6", "BARPA", "AWRA" )#CABLE", "AWRA")!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
datasets <- c("BARPA", "AWRA" )#CABLE", "AWRA")!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
    
    print(paste0("Plotting variable: ", vars[v], ", metric ", m, "/", length(metrics)))
    
    ### Set up figure ###
    png(paste0(outdir, "/FigureX", "_ensemble_changes_in_", metrics[m], "_",
               percentile, "_", scale, "_", vars[v], ".png"),
        height=7, width=8.3, units="in", res=400)
    
    
    par(mai=c(0.1, 0.1, 0.2, 0.1))
    par(omi=c(0.1, 0.3, 0.4, 0.1))
    
    # layout(matrix(c(1:3, 4, 4, 4, 5:7, 8, 8, 8), nrow=4, byrow=TRUE), heights=c(1, 0.3, 1, 0.3,
    #                                                                             1, 0.3, 1, 0.3,
    #                                                                             1, 0.3, 1, 0.3))
    # 
    
    par(mfrow=c(3, 3))
    par(bty="n")
    
    
    #####################
    ### Read datasets ###
    #####################

    #Loop through datasets    
    for (d in 1:length(datasets)) {
      
      
      plot_data_hist <- brick()
      plot_data_fut <- lapply(envelopes, function(x) brick())

      
      #Find models (use historical to find them, only want to include models that have historical data)
      #should be all of them but just in case
      #For AWRA/CABLE this finds the bc methods, doesn't matter for processing
      models <- list.files(paste0(path, "/Mean_drought_metrics/", datasets[d], 
                                       "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                                       "/", metrics[m], "/historical/"))
      
      
      for (mod in 1:length(models)) {
        
        mod_path <-  paste0(path, "/Mean_drought_metrics/", datasets[d], 
                            "/scale_", scale, "/", vars_all[[vars[v]]][[datasets[d]]], 
                            "/", metrics[m])
        


        ##################
        ### AWRA/CABLE ###
        ##################
        
        
        #Have bc methods in addition to GCMs
        if (datasets[d] %in% c("AWRA", "CABLE")) {
          
          #models correspond to BC methods here. Find GCMs to loop through
          gcms <- list.files(paste0(mod_path, "/historical/", models[mod]))
            
            
          #Initialise
          all_hist <- brick()
          all_fut  <- lapply(envelopes, function(x) brick())
          
          #Loop through ensembles
          for (gc in 1:length(gcms)) {

            #Find historical file
            hist_file <- list.files(paste0(mod_path, "/historical/", models[mod], "/",
                                           gcms[gc]), pattern=".nc", full.names=TRUE)
            
            #Should only find one file, check
            if (length(hist_file) != 1) stop("wrong CABLE/AWRA historical file")
            
            #Read data
            all_hist <- addLayer(all_hist, raster(hist_file))
            
            
            ### Future data ##
            
            for(gw in 1:length(envelopes)) {
              
              #1 degree of warming
              files_deg <- list.files(paste0(mod_path, "/", envelopes[gw], "deg/", models[mod], "/", 
                                             gcms[gc]),  recursive=TRUE, pattern=".nc", full.names=TRUE)
              
              #Check that found files, some GW levels don't have any
              if (length(files_deg) > 0) {
                
                #Read future data and calculate difference to historical
                data_deg <- brick(lapply(files_deg, raster)) - all_hist[[ens]]
                
                all_fut[[gw]] <- addLayer(all_fut[[gw]], calc(data_deg, median))
                
              }
              
            }

          } #GCMs            
            
            
          #Calculate median of ensemble members and resample to a common resolution
          ens_median_hist <- calc(all_hist, median)
          ens_median_fut  <- lapply(all_fut, function(x) if(!all(is.na(values(x)))) calc(x, median))
          
          
          #Collate
          plot_data_hist <- addLayer(plot_data_hist, ens_median_hist)
          
          plot_data_fut <- lapply(envelopes, function(x) if(!is.null(ens_median_fut[[x]]))
                                                                        addLayer(plot_data_fut[[x]],
                                                                        ens_median_fut[[x]])
                                                                        else plot_data_fut[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
          
          
        ##################
        ### CMIP6 data ###
        ##################
          
        } else if (datasets[d] %in% c("CMIP6", "BARPA")) {
          
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
          ens_median_fut  <- lapply(all_fut, function(x) if(!all(is.na(values(x)))) resample(calc(x, median), res_raster))
                
          
          #Collate
          plot_data_hist <- addLayer(plot_data_hist, ens_median_hist)

          plot_data_fut <- lapply(envelopes, function(x) if(!is.null(ens_median_fut[[x]]))
                                                                           addLayer(plot_data_fut[[x]],
                                                                                 ens_median_fut[[x]])
                                                 else plot_data_fut[[x]]) #this else just returns the same data, i.e. does nothing. Otherwise retursn a NULL and messes things up 
          
          #save model names
          
          
        } #CMIP6/BARPA
        
        
        
      } #models
      
      
      
      
      
      
      
      
      
      
      ################
      ### Plotting ###
      ################
      
      
      ### Historical median ###
      
      hist_median <- calc(plot_data_hist, median)
      
      #Plot
      len <- length(lims_hist[[vars[v]]][[metrics[m]]])
      image(hist_median, breaks=lims_hist[[vars[v]]][[metrics[m]]], 
            col=c(cols_hist(len-1)),
            axes=FALSE, ann=FALSE, asp=1)
      
      #Australia outline
      map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
      
      
      #time period label
      if (d == 1) mtext(side=3, "Historical mean")
        
      #dataset label
      mtext(side=2, datasets[d])
      
      
      
      
      ### Future difference ###
      
      for (env in 1:length(envelopes)) {
        
        #calculate future difference by model and then take median
        fut_median <- calc(plot_data_fut[[env]], fun=median)
        
        
        #Plot
        len <- length(lims_diff[[metrics[m]]])
        image(fut_median, breaks=lims_diff[[metrics[m]]], 
              col=c(cols_diff(len-1)),
              axes=FALSE, ann=FALSE, asp=1)
        
        #Australia outline
        map(region="Australia", add=TRUE, lwd=0.7) #border="grey50"
        
        
        #time period label
        if (d == 1) mtext(side=3, paste0(env, " degree"))
        
        
        
        
        
      }
        
      
      
      
        
    } #datasets
    
   
     

    
    # #Empty plot
    #  
    # plot(1, type="n", bty="n", yaxt="n", xaxt="n")
    # 
    # #Legend
    # len1 <- length(lims_diff[[metrics[m]]])-1
    # add_raster_legend2(cols=cols_diff(len1), limits=lims_diff[[metrics[m]]][2:len1],
    #                    main_title=unit[m], plot_loc=c(0.3,0.7,0.63, 0.77), 
    #                    title.cex=1, spt.cex=1, clip=TRUE, ysp_title_old=FALSE)
    # 
    
  
    dev.off ()
    
  } #metrics
  
} #variables






