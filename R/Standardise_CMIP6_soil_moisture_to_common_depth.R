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

#IPSL files don't include soil layer depths, can't use it
#remove file
files <- files[-which(grepl("mrsol//IPSL-", files))]


### Find maximum depth across models ###

#Open file handles
nc <- lapply(files, function(x) nc_open(x)) #tryCatch(nc_open(x), error=function(e) NA))

#Extract depth information (one model calls the bands by a different name, duh)
#depths <- lapply(nc, function(x) tryCatch(ncvar_get(x, "depth_bnds"), error=function(e) ncvar_get(x, "sdepth_bnds")))


#Writing as a loop so can deal with errors better
depths <- list()
for (k in 1:length(nc)) {
  
  depths[[k]] <- tryCatch(ncvar_get(nc[[k]], "depth_bnds"), error=function(e) tryCatch(ncvar_get(nc[[k]], "sdepth_bnds"),
                                                                                     error=function(e) ncvar_get(nc[[k]], "depth")))
}




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
      
      
     
      #Find file
      file <- list.files(paste(in_path, experiments[e], "mrsol", models[k], ensemble_members[ens], sep="/"),
                         full.names=TRUE, pattern="_Aus.nc")
      
      #Something wrong with ICON lat/lon bounds so it can't be processed by the main CMIP6
      #processing code. Skipping
      if (models[k] == "ICON-ESM-LR" & length(file) == 0) next
      
      #Need to skip IPSL because of insufficient soil layer information (see above)
      if (grepl("IPSL-", models[k])) next
      
      #Read in a layer as brick for extent information and timing
      data_raster <- brick(file, layer=1, level=1, stopIfNotEqualSpaced=FALSE)
      
      #Dataset years
      start_yr <- substr(names(data_raster)[1], 2,5)
      end_yr   <- substr(names(data_raster)[nlayers(data_raster)], 2,5)
      
      
      #Create output path
      out_path <- paste0(in_path, "/", experiments[e], "/mrsol_std_", max_depth, "m/", 
                         models[k], "/", ensemble_members[ens])
      dir.create(out_path, recursive=TRUE)
      
      out_file <- paste0(out_path, "/Monthly_standardised_mrsol_", max_depth, "m_", models[k], 
                         "_", ensemble_members[ens], "_", start_yr, "_", end_yr, "_Aus.nc")
      
      if (!file.exists(out_file)) {
        
       
        
        #Load all data
        nc <- nc_open(file)
        
        data <- ncvar_get(nc, "mrsol")
        
        #Get calendar and time units (need this for drought metric calculation)
        calendar <- ncatt_get(nc, 'time')$calendar
        
        tunits=ncatt_get(nc, 'time')$units
        
        time <- ncvar_get(nc, "time")
        
        
        
        #Initialise standardised data
        std_data <- array(NA, dim=c(dim(data)[1], dim(data)[2], dim(data)[4]))
        
        
        #Loop through years
        for(y in 1:dim(std_data)[3]){
          
          
          
          #Extract depth information
          depths <- tryCatch(ncvar_get(nc, "depth_bnds"), error=function(e) tryCatch(ncvar_get(nc, "sdepth_bnds"),
                                                                                     error=function(e) ncvar_get(nc, "depth")))
          
          #GFDL-ESM4 reports depths as a single vector
          #rather than depth bounds, need to convert depths into a
          #matrix (row x column ) to match other models
          if(is.na(ncol(depths))) {
            
            total_no <- length(depths)
            depths <- matrix(data=c(depths[1], rep(depths[2:(total_no-1)], each=2), depths[total_no]),
                             nrow=2, ncol=total_no-1)
          }
          
          
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
        
        # #output file name (create a temp one and a final one)
        # out_file_temp <- paste0(out_path, "/Monthly_standardised_mrsol_temp.nc")
        # 
        
        #Write output file
        
        
        xd = nc$dim[[3]]
        yd = nc$dim[[4]]
        
        # Define time dimension:
        td = nc$dim[[1]]
        
        
        sm_var <- ncdf4::ncvar_def(name=paste0("mrsol_std_", max_depth, "m"),
                                   units="kg m-2",
                                   dim=list(xd,yd,td),
                                   #missval=Nc_MissingVal,
                                   longname=paste("standardised mrsol", max_depth, "m"))
        
        # Create
        ncid <- ncdf4::nc_create(out_file, vars=sm_var)
        
        ncdf4::ncvar_put(nc=ncid,
                         varid=sm_var,
                         vals=std_data)
        
        
        # Close netcdf file:
        ncdf4::nc_close(ncid)
        
        
        
        # 
        # writeRaster(out_brick, out_file_temp, format="CDF", overwrite=TRUE, varname=paste0("mrsol_std_", max_depth, "m"), 
        #             longname=paste("standardised mrsol", max_depth, "m"), varunit="kg m-2",
        #             xname="longitude", yname="latitude", zname="time", zunit=tunits)
        # 
        
        rm(data_raster)
        rm(data)
        rm(std_data)
        
        # #Set calendar using CDO. Need this for drought calculation code so it can
        # #set time vector correctly
        # system(paste0("cdo setcalendar,", calendar, " -settime,", time, " -settunits,'", tunits, "' ", out_file_temp, " ", out_file))
        #   
        #file.remove(out_file_temp)
        
        
      }
      
    } #ensemble
  } #model
} #experiment













