# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:34:05 2016

@author: annaukkola

"""


from netCDF4 import Dataset,num2date # to work with NetCDF files
import numpy as np
import glob
import sys 
import os
import datetime
import xarray as xr


#fix to /g/data/w97/amu561/CMIP6_runoff_drought/CMIP6_data/Processed_CMIP6_data/ssp126/mrro/NorESM2-LM/r1i1p1f1
# file=mrro_NorESM2-LM_ssp126_r1i1p1f1_mm_month_2015_2100_ssp126_regrid_setgrid.n
# #Missing first January
# file=`ls *setgrid.nc`
# cp $file original_data_missing_tstep.nc
# cdo seltimestep,1 $file first_month.nc
# cdo shifttime,-1mon first_month.nc shifted.nc
# cdo mergetime shifted.nc $file fixed.nc
# mv fixed.nc $file
# rm shifted.nc first_month.nc

#All NorESM-LM and NorEMS-MM lai and mrro files were fixed
#Also /g/data/w35/amu561/CMIP6_runoff_drought//CMIP6_data/Processed_CMIP6_data/ssp585/mrro/TaiESM1/r1i1p1f1/
#mrro_TaiESM1_ssp585_r1i1p1f1_mm_month_2015_2100_ssp585_regrid_setgrid.nc
#and equivalent LAI file

### Set paths ###

root_path = "/g/data/w97/amu561/CABLE_AWRA_comparison/"
lib_path  = root_path + '/scripts/drought_scripts_AWRA/functions' # '/drought_metric_scripts/Drought_metrics/functions'

scratch_path = '/scratch/w97/amu561/'

sys.path.append(os.path.abspath(lib_path))
from drought_metrics import *

   
   
### Set variable ###

#Rainfall, total runoff, total soil moisture
#Using total soil moisture rather than standardising as JULES soil depth is
#3m which is close to the CMIP6 standardised depth
var_name=['pr', 'mrro', 'mrso']

var_path=var_name #This one is only used to create path/file names

#BoM data location
in_path="/g/data/ia39/australian-climate-service/release/CORDEX-CMIP6/output/AUS-15/BOM/"


#######################
### Set experiments ###
#######################

experiment=['historical', 'ssp126', 'ssp370'] 


#################
### Set years ###
#################

#Set baseline
baseline = [1970, 2005]


#Select if want indices returned as a time series or shorter vector
#If True, will return drought indices as time series (with NA for non-drought months)
#useful for calculating trends etc. in metrics
#If False, collapses indices into a short vector 
#reduces data size and useful if looking at the mean of drought indices
return_all_tsteps=True

############################
#### Set drought limits ####
############################

#Set percentile for drought threshold
perc=15

#Set scale if want to use running means to determine drought (like SPI scale)
#Use 1 if don't want to use this
scale=[3, 12]

#Use threshold determined separately for each month?
#If set to false, determines one threshold from all data.
#Set to false if using annual data
monthly=True


### Create temporary directory ###
temp_dir = str(scratch_path + "/monthly_sums_BARPA/")
os.system("mkdir -p " + temp_dir)


##################
### Load files ###
##################

#Loop through scales
for s in range(len(scale)):

    #List all model names
    models = os.listdir(in_path)
    
    #Remove ERA5 as only a historical evaluation dataset
    models.remove("ECMWF-ERA5")
 
    #Reference data for calculating threshold
    for m in range(len(models)):

        #Loop through experiments
        for k in range(len(experiment)):

            ensemble = os.listdir(in_path +  "/" + models[m] + "/" + experiment[k])

            #Loop through ensemble members (only one per model?)
            for e in range(len(ensemble)):
                
                #Loop through variables
                for v in range(len(var_name)):
        
                    #Progress
                    print("#--------------- Variable " + str(v) + "/" + str(len(var_name)))


                    #Use observations/historical data for calculating threshold?
                    #Uses this file to calculate baseline for drought metrics if true
                    #(currently set to use historical simulation, see "obs_file" below)
                    if experiment[k] == "historical":
                        obs_ref = False
                    else:
                        obs_ref = True
                        obs_var = var_name   #ET_mean, ET_min, ET_max


                    #Print progress
                    print('Variable: ' + var_name[v] + ', experiment: ' + experiment[k] + 
                          '(' + str(k+1) + '/' + str(len(experiment)) + '), model: ' + 
                          models[m] + '(' + str(m+1) + "/" + str(len(models)))
                
                                             
                    
                    ### Find BARPA files ###
                    
                    #List files
                    files = glob.glob(in_path +  "/" + models[m] + "/" + experiment[k] + "/" +
                               ensemble[e] + "/BOM-BARPA-R/v1/mon/" + var_name[v] + "/*.nc")
                    
                      
                    #If using obs as reference
                    if obs_ref:
                        obs_files = glob.glob(in_path +  "/" + models[m] + "/historical/" +
                                   ensemble[e] + "/BOM-BARPA-R/v1/mon/" + var_name[v] + "/*.nc")
                        

                        


                    ########################
                    ### Pre-process data ###
                    ########################
                    
                    
                    #Temporary file for saving monthly totals
                    temp_file = str(temp_dir + "/" + experiment[k] + "_" + 
                               models[m] + "_" + ensemble[e] + "_" + var_name[v] + ".nc")


                    #Calculate monthly sums if file doesn't already exist       
                    if not os.path.isfile(temp_file):

                            
                        #Model data
                        ds = xr.open_mfdataset(files)
                        
                        data = ds[var_name[v]]                        
                        
                        
                        #For runoff and precip, multiple by the number of days per month and
                        #seconds per day to go from mm/s to mm/month
                        if var_name[v] == "mrro" or var_name[v] == "pr":
                            data = data * data.time.dt.daysinmonth * 86400.0
                        
                        
                        if var_name[v] == "pr":
                            
                            #Need to read in runoff file to mask out ocean grid cells
                            mrro_file = files = glob.glob(in_path +  "/" + models[m] + "/" + experiment[k] + "/" +
                                       ensemble[e] + "/BOM-BARPA-R/v1/mon/mrro/*.nc")
                            
                            #Read runoff data
                            mrro_ds = xr.open_mfdataset(mrro_file)
                            
                            mrro = mrro_ds.mrro                       
                            
                            #Mask where runoff is NA
                            data = xr.where(np.isnan(mrro),np.nan, data)
                            
                        
                        #Crop to exclude regions outside Australia
                        #Use CABLE/AWRA dimensions                        
                        data = data.sel(lat=slice(-44,-10), lon=slice(112,154))

                        
                        #Rename data variable
                        data.name = var_name[v]
                        

                        #Write to file
                        data.to_netcdf(temp_file, format='NETCDF4', 
                                       encoding={var_name[v]:{
                                                 'shuffle':True,
                                                 'chunksizes':[12, 200, 40],
                                                 'zlib':True,
                                                 'complevel':5}
                                                 })


                    #If using obs as reference
                    if obs_ref:
                        
                        #Temporary file for saving monthly totals
                        temp_ref_file = str(temp_dir + "/historical_" + models[m] + 
                                            "_" + ensemble[e] + "_" + var_name[v] + ".nc")


                        #Calculate monthly sums if file doesn't already exist       
                        if not os.path.isfile(temp_file):


                            #Model data
                            ds_ref = xr.open_mfdataset(obs_files)
                            
                            data_ref = ds_ref[var_name[v]]
                            
                            #For runoff and precip, multiple by the number of days per month and
                            #seconds per day to go from mm/s to mm/month
                            if var_name[v] == "mrro" or var_name[v] == "pr":
                                data_ref = data_ref * data_ref.time.dt.daysinmonth * 86400.0
                                                    

                            if var_name[v] == "pr":
                                data_ref = xr.where(np.isnan(mrro),np.nan, data_ref)
                                
                            
                            #Crop to exclude regions outside Australia
                            #Use CABLE/AWRA dimensions                        
                            data_ref = data_ref.sel(lat=slice(-44,-10), lon=slice(112,154))


                            #Rename data variable
                            data_ref.name = var_name[v]
                            

                            #Write to file
                            data_ref.to_netcdf(temp_ref_file, format='NETCDF4', 
                                               encoding={var_name[v]:{
                                                     'shuffle':True,
                                                     'chunksizes':[12, 200, 40],
                                                     'zlib':True,
                                                     'complevel':5}
                                                     })


                    #################
                    ### Load data ###
                    #################

                                        
                    fh       = Dataset(temp_file, mode='r')
                    all_data = fh.variables[var_name[v]][:]

                    data     = all_data.data 
                    fh_time  = fh.variables["time"]
                    mask     = all_data.mask
                    
                    #Get lon and lat 
                    lat = fh.variables['lat'][:] #northing
                    lon = fh.variables['lon'][:] #easting

                                     
                    #Get dataset dates
                    fh_dates = num2date(fh_time[:], fh_time.units, calendar=fh_time.calendar)
                    fh_years = np.array([y.year for y in fh_dates])


                    #fh.close()
                        
                    # #Mask
                    # #If mask one value, make it an array that has same dimensions
                    # #as data (I think this happens because oceans haven't been
                    # #masked out)
                    # #Some files don't have mask at all, hmmm
                    # try:
                    #     mask = all_data.mask                    
                    # except: 
                    #     mask = np.zeros(data.shape, dtype='bool')
                    # 
                    # if mask.shape == ():
                    #     mask = np.zeros(data.shape, dtype='bool')
                                        
                                    
                    miss_val = -99999.0
                    data[mask==True] = miss_val

                    miss_val = -99999.0
                    data[mask==True] = miss_val
                
                    #Read reference data used to calculate threshold
                    if obs_ref:
                        obs_fh      = Dataset(temp_ref_file, mode='r')            
                        control_ref = fh.variables[var_name[v]][:].data
                        #control_ref     = all_data.data 
                        obs_time   = fh.variables["time"]
                        
                        #Get lon and lat 
                        lat_ctrl = fh.variables['lat'][:] #northing
                        lon_ctrl = fh.variables['lon'][:] #easting

                        print("Using an OBS/model ref to calculate baseline")

                        #Get dataset dates
                        obs_dates = num2date(obs_time[:], obs_time.units,
                                             calendar=obs_time.calendar)

                        obs_years = np.array([y.year for y in obs_dates])


                    else:
                        control_ref = data        
                        
                        
                    
                    #Not sure why but python reads some files upside down
                    #Flip latitude if reads map upside down so model matches SPI data
                    if obs_ref:
                        #If model data upside down
                        if (lat[0] < 0 and lat_ctrl[0] > 0):
                            print("Flipping MODEL data, experiment: ", experiment[k], 
                                  " model:", models[m])
                            data = data[:,::-1,:]
                            #replace lat with lat_ctrl (otherwise written to file the wrong way round)
                            lat=lat_ctrl
                        #If REF data upside down
                        elif (lat[0] > 0 and lat_ctrl[0] < 0):
                            print("Flipping OBS ref data, experiment: ", experiment[k], 
                                  " model:", models[m])
                            control_ref = control_ref[:,::-1,:]


                    
                    ###################################################
                    ### Create output file name and check if exists ###
                    ###################################################
                    
                    #Creating this here so can check if it already exists,
                    #and skip to speed up processing
                            
                    #Creat output path
                    out_path = (root_path + '/Drought_metrics/BARPA/' + experiment[k] + "/" + 
                                var_path[v] ) # + '/Obs_' + str(obs_ref) )

                    #Add percentile
                    out_path = (out_path + '/Perc_' + str(perc) + 
                                '/Baseline_' + str(baseline[0]) + "_" + str(baseline[1]) +
                                "/Scale_" + str(scale[s]) + "/" + models[m] + "/" + ensemble[e])
                                
                    
                    #Create output directory if doesn't exist
                    if not os.path.exists(out_path):    
                        os.makedirs(out_path)
                    
                    #Create output file name
                    out_file = (out_path + '/' + models[m] + "_" + ensemble[e] + 
                                '_drought_metrics_perc_' + str(perc))

                    out_file = (out_file + "_" + experiment[k] + '_' + str(fh_years[0]) + 
                                '_' + str(fh_years[-1]) + '.nc')
                            
                    
                    #Check if exists and Skip
                    if os.path.isfile(out_file):
                        print("Skipping " + experiment[k] + ", " + models[m] + ", " + 
                              ensemble[e] + ", already exists")
                        continue
                    

                    #############################
                    ### Find baseline indices ###
                    #############################
                    
                    #Get dates
                    ref_years = fh_years
                    if obs_ref:
                        ref_years = obs_years
                        
                    #Find indices corresponding to start of first year,
                    #and end of end year, defined by baseline
                    
                    #Create indices for baseline period
                    subset = range(np.where(ref_years == baseline[0])[0][0],
                                   np.where(ref_years == baseline[1])[0][-1] + 1) #Stupid python indexing
                                    


                    ################################
                    ### Initialise output arrays ###
                    ################################
                    
                    #Expect to have approx. percentile amount of months as drought (e.g. 10% of data
                    #when percentile is set to 10. To be on the safe side, determine array sizes
                    #as double that size) if not writing out as a full time series
                    
                    if return_all_tsteps:
                        save_len = len(data)
                    else:
                        save_len = int(len(data)*(perc_onset/100)*2)
                
                    duration          = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
                    rel_intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
                    rel_intensity_mon = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan

                    timing            = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    

                    if monthly:
                        threshold         = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan

                    else:
                        threshold         = np.zeros((len(lat), len(lon))) + miss_val # * np.nan
                
                
                    #########################
                    ### Calculate metrics ###
                    #########################


                    #Loop through grid cells
                    for i in range(len(lat)):
                                
                        for j in range(len(lon)):
                        
                             #Calculate metrics if cell not missing
                             if (any(~mask[:,i,j]) and 
                             len(np.unique(data[:,i,j][~np.isnan(data[:,i,j])])) > 1):  
                             
                                 #Calculate metrics
                                 metric = drought_metrics(mod_vec=data[:,i,j], lib_path=lib_path,
                                                          perc=perc, miss_val=miss_val,
                                                          monthly=monthly, obs_vec=control_ref[:,i,j],
                                                          return_all_tsteps=return_all_tsteps, scale=scale[s],
                                                          add_metrics=(['rel_intensity', 'threshold', 
                                                          'rel_intensity_monthly','timing']),
                                                          subset=subset)
                            
        
                                 ### Write metrics to variables ###
                                 duration[range(np.size(metric['duration'])),i,j] = metric['duration']  #total drought duration (months)

                                 rel_intensity[range(np.size(metric['rel_intensity'])),i,j] = metric['rel_intensity'] #average magnitude

                                 rel_intensity_mon[range(np.size(metric['rel_intensity_monthly'])),i,j] = metric['rel_intensity_monthly'] #average magnitude
                            
                                 timing[range(np.size(metric['timing'])),i,j] = metric['timing']    #drought timing (month index)

                                 if monthly:
                                     threshold[:,i,j] = metric['threshold'][0:12]    #drought timing (month index)
                                 else:
                                     threshold[i,j]   = metric['threshold']
                        
                                
                    ##############################
                    ### Write result to NetCDF ###
                    ##############################
                                
                    # Open a new netCDF file for writing
                    ncfile = Dataset(out_file,'w', format="NETCDF4_CLASSIC") 
                    
                    # Create the output data
                    # Create the x, y and time dimensions
                    ncfile.createDimension('lat', lat.shape[0])
                    ncfile.createDimension('lon', lon.shape[0])
                    ncfile.createDimension('time', save_len)
                    
                    
                    if monthly:
                        ncfile.createDimension('month', 12)

                    # Create dimension variables
                    longitude = ncfile.createVariable("lon",  'f8', ('lon',))
                    latitude  = ncfile.createVariable("lat",  'f8', ('lat',))
                    time      = ncfile.createVariable("time", 'f8', ('time',))

                    if monthly:
                        month = ncfile.createVariable("month", 'i4', ('month',))

                    #Create data variables
                    data_dur  = ncfile.createVariable('duration', 'f8',('time','lat','lon'), fill_value=miss_val)
                    data_mag  = ncfile.createVariable('rel_intensity','f8',('time','lat','lon'), fill_value=miss_val)
                    data_rel  = ncfile.createVariable('rel_intensity_by_month','f8',('time','lat','lon'), fill_value=miss_val)
                    data_tim  = ncfile.createVariable('timing', 'i4',('time','lat','lon'), fill_value=miss_val)

                    #Create data variable for threshold
                    if monthly:
                        data_thr = ncfile.createVariable('threshold', 'f8',('month','lat','lon'), fill_value=miss_val)
                    else:
                        data_thr = ncfile.createVariable('threshold', 'f8',('lat','lon'), fill_value=miss_val)


                    #Set variable attributes
                    longitude.units = 'degrees_east'
                    latitude.units  = 'degrees_north'
                    time.units      = fh_time.units

                    #Calendar
                    time.calendar   = fh_time.calendar

                    data_dur.long_name = 'drought event duration (no. months)'
                    data_mag.long_name = 'drought event relative intensity (%)'
                    data_rel.long_name = 'drought month relative intensity (%)'
                    data_tim.long_name = 'drought event timing (binary drought/non-drought index)'
                    data_thr.long_name = 'drought threshold (mm)'

                    if monthly:
                        month[:] = range(1,12+1)

                    # Write data to dimension variables
                    longitude[:]=lon
                    latitude[:] =lat

                    #If saving all time steps
                    if return_all_tsteps:
                        time[:] = fh_time[:]
                    else:
                        time[:] = range(1, save_len+1)
                            
                    if monthly:
                        month[:] = range(1,12+1)

                    #Write data to data variables
                    data_dur[:,:,:] = duration    
                    data_mag[:,:,:] = rel_intensity
                    data_rel[:,:,:] = rel_intensity_mon
                    data_tim[:,:,:] = timing

                    if monthly:    
                        data_thr[:,:,:] = threshold
                    else:
                        data_thr[:,:]   = threshold

                    # Close the file
                    ncfile.close()


                    #Finally compress file
                    #Don't use overwrite option as fails for some files
                    os.system("nccompress " + out_file) 
                    os.system("mv " + out_path + "/tmp.nc_compress/*.nc" + " " + out_file)
                    os.system("rm -r " + out_path + "/tmp.nc_compress/")
