library(ncdf4)
library(raster)
library(lubridate)

# initial garbage collection
rm( list=ls(all=TRUE) )

###--- Function ---###
in_interval <- function(x, interval){ 
  stopifnot(length(interval) == 2L) 
  interval[1] < x & x <= interval[2] 
} 


#################
### Set paths ###
#################

path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"

in_path <- paste0(path, "/CMIP6_data/Processed_CMIP6_data/")


experiments <- c("historical", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585")

#Find models
models <- list.files(paste0(in_path, "/historical/mrsol"))



#Find mrso files
files <- list.files(paste0(in_path, "/historical/mrsol/"), pattern="_Aus.nc", 
                    full.names=TRUE, recursive=TRUE)


### Find maximum depth across models ###

#Open file handles
nc <- lapply(files, function(x) nc_open(x)) #tryCatch(nc_open(x), error=function(e) NA))

#Extract depth information (one model calls the bands by a different name, duh)
depths <- lapply(nc, function(x) tryCatch(ncvar_get(x, "depth_bnds"), error=function(e) ncvar_get(x, "sdepth_bnds")))



#Find common max depth across datasets
max_depth <- min(sapply(depths, max), na.rm=TRUE)


lapply(nc, nc_close)




### Calculate standardised soil moisture ###


#Loop through experiments
for (e in 1:length(experiments)) {
  
  
  models <- list.files(paste(in_path, experiments[e], "mrsol", sep="/"))
  
  if(length(models) == 0) next
  
  #Loop through models
  for(k in 1:length(models)){
    
    ensemble_members <- list.files(paste(in_path, experiments[e], "mrsol", models[k], sep="/"))
    
    
    #Loop through ensemble members
    for (ens in 1:length(ensemble_members)) {
      
      
      #Create output path
      out_path <- paste0(in_path, "/", experiments[e], "/mrsol_std_", max_depth, "m/", 
                         models[k], "/", ensemble_members[ens])
      dir.create(out_path, recursive=TRUE)
      
    
      #Find file
      file <- list.files(paste(in_path, experiments[e], "mrsol", models[k], ensemble_members[ens], sep="/"),
                         full.names=TRUE, pattern="_Aus.nc")
      
      
      #Read in a layer as brick for extent information and timing
      data_raster <- brick(file, layer=1, level=1, stopIfNotEqualSpaced=FALSE)
      
      
      #Load all data
      nc <- nc_open(file)
      
      data <- ncvar_get(nc, "mrsol")
      
      
      #Initialise standardised data
      std_data <- array(NA, dim=c(dim(data)[1], dim(data)[2], dim(data)[4]))
      
      
      #Loop through years
      for(y in 1:dim(std_data)[3]){
        
        
        
        #Extract depth information
        depths <- tryCatch(ncvar_get(nc, "depth_bnds"), error=function(e) ncvar_get(nc, "sdepth_bnds"))
        
        
        #Find which layer encompasses max_depth
        max_layer <- which(sapply(1:ncol(depths), function(y)  in_interval(max_depth, depths[,y])))
        
        
        #Sum all layers above max layer
        std_data[,,y] <- apply(data[,,1:(max_layer-1), y], MARGIN=c(1,2), sum)

        
        #Depth of maximum layer
        depth_maxl <- depths[2,max_layer]-depths[1,max_layer]
        
        #Difference of max depth and lower bound of max layer
        depth_std <- max_depth - depths[1,max_layer]
        
        #Calculate weight (proportion of max layer covered by max depth)
        weight <- depth_std / depth_maxl
        
        
        #Add weighted max layer to standardised data      
        std_data[,,y] <- std_data[,,y] + data[,,max_layer, y] * weight
        
      }
      
      nc_close(nc)
      
      
      #Convert to raster brick
      
      out_brick <- brick(std_data)
      
      #Flip if not matching raster
      if(nrow(out_brick) != nrow(data_raster)){
        out_brick <- flip(t(out_brick), direction='y')
      }
      
      #Set extent and layer names
      
      
      
      extent(out_brick) <- extent(data_raster)
      names(out_brick)  <- names(data_raster)
      
      ### Write output ###
      
      start_yr <- substr(names(out_brick)[1], 2,5)
      end_yr   <- substr(names(out_brick)[nlayers(out_brick)], 2,5)
      
      out_file <- paste(out_path, "/Monthly_standardised_mrsol_", max_depth, "m_", models[k], 
                        "_", ensemble_members[ens], "_", start_yr, "_", end_yr, "_Aus.nc", sep="")
      
      
      #Write output file
      writeRaster(out_brick, out_file, format="CDF", overwrite=TRUE, varname=paste0("mrsol_std_", max_depth, "m"), 
                  longname=paste("standardised mrsol", max_depth, "m"), varunit="kg m-2",
                  xname="longitude", yname="latitude", zname="time", zunit=paste("months since Jan", start_yr))
      
      
      rm(data_raster)
      rm(data)
      rm(std_data)
      
        
    } #ensemble
  } #model
} #experiment













