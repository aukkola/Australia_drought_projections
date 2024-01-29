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


#Set scale if want to use running means to determine drought (like SPI scale)
#Use 1 if don't want to use this
scale=[3]



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


        #Check that all experiments exist. Some models only have a subset
        #available, skip these
        if len(os.listdir(in_path +  "/" + models[m] )) != len(experiment):
            continue
            
        #These models are not yet available even though folder exists, skip
        if models[m] == "CMCC-CMCC-ESM2" or models[m] == "NCAR-CESM2":
            continue

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
                    out_path = (root_path + '/BARPA_data/' + experiment[k] + "/" + 
                                var_path[v] +"/" + models[m] + "/" + ensemble[e])
                                #"/Scale_" + str(scale[s]) + "/" + models[m] + "/" + ensemble[e])
    
                    os.system("mkdir -p " + out_path)


                    #Temporary file for saving monthly totals
                    temp_file = str(out_path + "/" + experiment[k] + "_" + 
                               models[m] + "_" + ensemble[e] + "_" + var_name[v] + 
                               "_running_mean_" + str(scale[s]) + "months.nc")


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

                        #Calculate rolling mean time series
                        #data = data.rolling(time=scale[s], center=True).mean()
                        
                        
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


                    # #If using obs as reference
                    # if obs_ref:
                    # 
                    #     #Temporary file for saving monthly totals
                    #     temp_ref_file = str(temp_dir + "/historical_" + models[m] + 
                    #                         "_" + ensemble[e] + "_" + var_name[v] + ".nc")
                    # 
                    # 
                    #     #Calculate monthly sums if file doesn't already exist       
                    #     if not os.path.isfile(temp_file):
                    # 
                    # 
                    #         #Model data
                    #         ds_ref = xr.open_mfdataset(obs_files)
                    # 
                    #         data_ref = ds_ref[var_name[v]]
                    # 
                    #         #For runoff and precip, multiple by the number of days per month and
                    #         #seconds per day to go from mm/s to mm/month
                    #         if var_name[v] == "mrro" or var_name[v] == "pr":
                    #             data_ref = data_ref * data_ref.time.dt.daysinmonth * 86400.0
                    # 
                    # 
                    #         if var_name[v] == "pr":
                    #             data_ref = xr.where(np.isnan(mrro),np.nan, data_ref)
                    # 
                    # 
                    #         #Crop to exclude regions outside Australia
                    #         #Use CABLE/AWRA dimensions                        
                    #         data_ref = data_ref.sel(lat=slice(-44,-10), lon=slice(112,154))
                    # 
                    #         #Calculate rolling mean time series
                    #         data_ref = data_ref.rolling(time=scale, center=True).mean()
                    # 
                    # 
                    #         #Rename data variable
                    #         data_ref.name = var_name[v]
                    # 
                    # 
                    #         #Write to file
                    #         data_ref.to_netcdf(temp_ref_file, format='NETCDF4', 
                    #                            encoding={var_name[v]:{
                    #                                  'shuffle':True,
                    #                                  'chunksizes':[12, 200, 40],
                    #                                  'zlib':True,
                    #                                  'complevel':5}
                    #                                  })
