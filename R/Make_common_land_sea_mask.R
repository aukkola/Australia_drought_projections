#Use AWRA input file to make a common 5km land-sea mask

library(raster)

#Set path
path <- "/g/data/w97/amu561/CABLE_AWRA_comparison/"


#Get AWRA file
awra_file <- paste0(path, "/..//Steven_CABLE_runs/AWRA_fields",
                     "/AWRA_fractional_vegetation_cover_monthly_climatology_1960_2005.nc")
                     
#Read data
data <- raster(awra_file)

#Set all non-NA values to 1
data[!is.na(data)] <- 1
data[is.na(data)]  <- 0


#Save data
outdir <- paste0(path, "/landsea/")
dir.create(outdir)

writeRaster(data, paste0(outdir, "/landsea_mask_AWRA.nc"), varname="landsea", overwrite=TRUE)



#!!!! NB !!!!!

#Need to fix the grid afterwards for this to work with CDO
#No idea why there is an issue in the first place

#Run this on the landsea file
#cdo griddes landsea_mask_AWRA.nc > grid.txt

#Then change the grid.txt so that it looks like this:

## 
## gridID 1
## 
# gridtype  = lonlat
# gridsize  = 572721
# xsize     = 841
# ysize     = 681
# xname     = longitude
# xlongname = "longitude"
# xunits    = "degrees_east"
# yname     = latitude
# ylongname = "latitude"
# yunits    = "degrees_north"
# xfirst    = 112
# xinc      = 0.05
# yfirst    = -10
# yinc      = -0.05

#Then set the new grid:
#cdo setgrid,grid.txt landsea_mask_AWRA.nc landsea_mask_AWRA_fixed.nc 







