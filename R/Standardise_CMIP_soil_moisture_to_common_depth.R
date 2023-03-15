library(ncdf4)
library(raster)


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

datasets <- c("CMIP6", "CMIP5")


for (d in datasets) {
  
  in_path <- 
  
  
  
  
  
  
  
  
}




experiment <- "historical"

#Find models
models <- list.files(paste(path, experiment, "mrlsl", sep="/"))

#Find mrlsl files
files <- lapply(models, function(x) list.files(paste(path, experiment, "mrlsl", x, sep="/"), 
                                               pattern="regrid.nc", full.names=TRUE))


### Find depth informaton ###

#Open file handles
nc <- lapply(files, function(x) nc_open(x)) #tryCatch(nc_open(x), error=function(e) NA))

#Extract depth information
depths <- lapply(nc, function(x) ncvar_get(x, "depth_bnds")) #tryCatch(ncvar_get(x, "depth_bnds"), error=function(e) NA))

#Find common max depth across datasets
max_depth <- min(sapply(depths, max), na.rm=TRUE)

#Find which layer encompasses max_depth
max_layer <- lapply(depths, function(x) if(all(!is.na(x))) which(sapply(1:ncol(x), function(y) 
  in_interval(max_depth, x[,y])))
  else NA)



### Calculate standardised soil moisture ###

for(k in 1:length(files)){
  
  
  #Create output path
  out_path <- paste(path, experiment, "std_mrso", models[k], sep="/")
  dir.create(out_path, recursive=TRUE)
  
  
  if(!is.na(depths[[k]])){
    
    #Read in a layer as brick for extent information and timing
    data_raster <- brick(files[[k]], layer=1, level=1, stopIfNotEqualSpaced=FALSE)
    
    
    #Load all data
    data <- ncvar_get(nc[[k]], "mrlsl")
    
    
    #Initialise standardised data
    std_data <- array(NA, dim=c(dim(data)[1], dim(data)[2], dim(data)[4]))
    
    
    #Loop through years
    for(y in 1:dim(std_data)[3]){
      
      #Sum all layers above max layer
      std_data[,,y] <- apply(data[,,1:(max_layer[[k]]-1), y], MARGIN=c(1,2), sum)
      
      #Depth of maximum layer
      depth_maxl <- depths[[k]][,max_layer[[k]]][2]-depths[[k]][,max_layer[[k]]][1]
      
      #Difference of max depth and lower bound of max layer
      depth_std <- max_depth - depths[[k]][,max_layer[[k]]][1]
      
      #Calculate weight (proportion of max layer covered by max depth)
      weight <- depth_std / depth_maxl
      
      
      #Add weighted max layer to standardised data      
      std_data[,,y] <- std_data[,,y] + data[,,max_layer[[k]], y] * weight
      
    }
    
    
    
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
    end_yr <- substr(names(out_brick)[nlayers(out_brick)], 2,5)
    
    out_file <- paste(out_path, "/Monthly_standardised_mrso_", max_depth, "m_", models[k], 
                      "_", start_yr, "_", end_yr, "_regrid.nc", sep="")
    
    
    
    writeRaster(out_brick, out_file, format="CDF", overwrite=TRUE, varname="std_mrso", 
                longname=paste("standardised mrlsl", max_depth, "m"), varunit="kg m-2",
                xname="longitude", yname="latitude", zname="time", zunit=paste("months since Jan", start_yr))
    
    
    rm(data_raster)
    rm(data)
    rm(std_data)
    
    
    
  }
}













