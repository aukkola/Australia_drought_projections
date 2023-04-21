              library(raster)
              library(ncdf4)


              ### Map of mean of all time slices ###
              
              files_regrid <- list.files(path="/g/data/w97/amu561/CABLE_AWRA_comparison/CMIP6_data/Processed_CMIP6_data//historical/pr/EC-Earth3-Veg-LR/", recursive=TRUE, 
                                         pattern="setgrid.nc", full.names=TRUE)    #regridded
                                         
              data_regrid <- lapply(files_regrid, brick, stopIfNotEqualSpaced=FALSE)


              pdf("/g/data/w97/amu561/CABLE_AWRA_comparison/CMIP6_data/Processed_CMIP6_data//Data_checks/Plots/historical/EC-Earth3-Veg-LR/pr//pr_historical_EC-Earth3-Veg-LR_r1i1p1f1_monthly_mean_regridded.pdf", 
                  height=5, width=8)
              par(mai=c(0.2,0.2,0.2,0.6))              
							par(omi=c(0.2,0.2,0.2,0.6))

              par(mfcol=c(ceiling(sqrt(length(data_regrid))), ceiling(sqrt(length(data_regrid)))))
              
              lapply(data_regrid, function(x) plot(mean(x), ylab="", xlab="", yaxt="n", xaxt="n"))
              dev.off()


          		### Global mean time series ###
              
          		files_mean <- list.files(path="/g/data/w97/amu561/CABLE_AWRA_comparison/CMIP6_data/Processed_CMIP6_data//Data_checks/historical/pr/Global_mean/EC-Earth3-Veg-LR/", recursive=TRUE, 
                                       pattern="EC-Earth3-Veg-LR_global_mean", full.names=TRUE)

          		nc_handles <- lapply(files_mean, nc_open)
          		nc_data    <- lapply(nc_handles, ncvar_get, varid="pr")

          		pdf("/g/data/w97/amu561/CABLE_AWRA_comparison/CMIP6_data/Processed_CMIP6_data//Data_checks/Plots/historical/EC-Earth3-Veg-LR/pr//pr_historical_EC-Earth3-Veg-LR_r1i1p1f1_global_mean_timeseries.pdf", 
                  height=13, width=40)
          		par(mai=c(0.2,0.2,0.2,0.6))
							par(omi=c(0.4,0.6,0.2,0.6))
          		par(mfcol=c(ceiling(sqrt(length(nc_data))), ceiling(sqrt(length(nc_data)))))
              
          		lapply(nc_data, function(x) plot(x, type="l", col="blue", 
										                           ylab=pr, xlab="Time step"))
          		dev.off()


