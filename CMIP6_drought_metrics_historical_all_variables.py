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


### Set paths ###

root_path = "/g/data/w97/amu561/CABLE_AWRA_comparison/"
lib_path  = root_path + '/Drought_metric_codes/functions' # '/drought_metric_scripts/Drought_metrics/functions'


sys.path.append(os.path.abspath(lib_path))
from drought_metrics import *

   
   
### Set variable ###

#Rainfall, total runoff, standardised soil moisture, total soil moisture
var_name=['pr' 'mrro' 'mrsol_std', 'mrso']

var_path=var_name #This one is only used to create path/file names


#######################
### Set experiments ###
#######################

experiment=['historical'] 


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



#Use observations/historical data for calculating threshold?
#Uses this file to calculate baseline for drought metrics if true
#(currently set to use historical simulation, see "obs_file" below)
obs_ref = False
obs_var = var_name   #ET_mean, ET_min, ET_max



##################
### Load files ###
##################

#Loop through scales
for s in range(len(scale)):


    #Loop through variables
    for v in range(len(var_name)):
        
        #Progress
        print("#--------------- Variable " + str(v) + "/" + str(len(var_name)))


        #Loop through experiments
        for k in range(len(experiment)):

            #List all model names
            models = os.listdir(root_path + '/CMIP6_data/Processed_CMIP6_data/' + 
                                experiment[k] + "/" + var_name[v] + "/")


            #Reference data for calculating threshold
            for m in range(len(models)):

                ensembles=os.listdir(root_path + '/CMIP6_data/Processed_CMIP6_data/' + 
                                     experiment[k] + "/" + var_name[v] + "/" + models[m])
                
                for e in range(len(ensembles)):
                
                
                    #Print progress
                    print('Experiment: ' + str(k+1) + '/' + str(len(experiment)) + ', model: ' + 
                          str(m+1) + "/" + str(len(models)) + ', ensemble:' + str(e+1) + '/' + 
                          str(len(ensembles)) + " (model name: " + models[m] + ")")
            
                      
                    #If using obs as reference (set to ET data, fix later if needed...)
                    if obs_ref:
                        obs_file = glob.glob(root_path + '/CMIP6_data/Processed_CMIP6_data/historical/' +
                                             var_name[v] + "/" + models[m] + "/" + ensembles[e] +
                                             "/*_Aus.nc")
                    
                
                    ### Find CMIP5 files ###
                    files = glob.glob(root_path + '/CMIP6_data/Processed_CMIP6_data/' + experiment[k] + 
                                      "/" + var_name[v] + "/" + models[m] + "/" + ensembles[e] + 
                                      "/*_Aus.nc")
                    
                    
                    #Skip EC-Earth3 runoff, missing years. Fix when available !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    #if var_name[v]=="mrro" and models[m]=="EC-Earth3" and ensembles[e]=="r1i1p1f1":
                    #    print("Skipping EC-Earth3 runoff, files corrupt !!!!!!!!!!!!!!!!!!!!!!!!!")
                    #    continue


                    #################
                    ### Load data ###
                    #################
                    
                    #Model data
                    fh = Dataset(files[0], mode='r')
                    all_data = fh.variables[var_name[v]][:] #[yr_ind]
                    data     = all_data #.data #####TEMPORARY whilst oceans not masked out !!!!!!!!!!!!!!
                    #mask     = all_data.mask
                    fh_time     = fh.variables["time"]

                    #Get lon and lat (name varies by CMIP5 model)
                    try:
                        lat = fh.variables['latitude'][:]
                        lon = fh.variables['longitude'][:]
                    except:
                        lat = fh.variables['lat'][:] #northing
                        lon = fh.variables['lon'][:] #easting

                 
                    #Get dataset dates
                    fh_dates = num2date(fh_time[:], fh_time.units, calendar=fh_time.calendar)
                    fh_years = np.array([y.year for y in fh_dates])


                    #fh.close()
                    
                    #Mask
                    #If mask one value, make it an array that has same dimensions
                    #as data (I think this happens because oceans haven't been
                    #masked out)
                    #Some files don't have mask at all, hmmm
                    try:
                        mask = all_data.mask                    
                    except: 
                        mask = np.zeros(data.shape, dtype='bool')

                    if mask.shape == ():
                        mask = np.zeros(data.shape, dtype='bool')
                                    
                                    
                    miss_val = -99999.0
                    data[mask==True] = miss_val

                    miss_val = -99999.0
                    data[mask==True] = miss_val
                
                    #Read reference data used to calculate threshold
                    if obs_ref:
                        obs_fh      = Dataset(obs_file[0], mode='r')
                        control_ref = obs_fh.variables[obs_var[v]][:] #.data ### TEMPORARY !!!!!!!!!!!!!!!!!!
                        obs_time    = obs_fh.variables["time"]
                    
                        #Get lon and lat (name varies by CMIP5 model)
                        try:
                            lat_ctrl = fh.variables['latitude'][:]
                            lon_ctrl = fh.variables['longitude'][:]
                        except:
                            lat_ctrl = fh.variables['lat'][:]
                            lon_ctrl = fh.variables['lon'][:]
                         

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
                    out_path = (root_path + '/Drought_metrics/CMIP6/' + experiment[k] + "/" + 
                                var_path[v] ) # + '/Obs_' + str(obs_ref) )

                    #Add percentile
                    out_path = (out_path + '/Perc_' + str(perc) + 
                                '/Baseline_' + str(baseline[0]) + "_" + str(baseline[1]) +
                                "/Scale_" + str(scale[s]) + "/" + models[m] + "/" + ensembles[e])
                                
                    
                    #Create output directory if doesn't exist
                    if not os.path.exists(out_path):    
                        os.makedirs(out_path)
                    
                    #Create output file name
                    out_file = (out_path + '/' + models[m] + "_" + ensembles[e] + 
                                '_drought_metrics_perc_' + str(perc))

                    out_file = (out_file + "_" + experiment[k] + '_' + str(fh_years[0]) + 
                                '_' + str(fh_years[-1]) + '.nc')
                            
                    
                    #Check if exists and Skip
                    if os.path.isfile(out_file):
                        print("Skipping " + experiment[k] + ", " + models[m] + ", " + 
                              ensembles[e] + ", already exists")
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

                        timing        = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    

                    if monthly:
                        threshold_onset       = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan
                        threshold_termination = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan

                    else:
                        threshold_onset       = np.zeros((len(lat), len(lon))) + miss_val # * np.nan
                        threshold_termination = np.zeros((len(lat), len(lon))) + miss_val # * np.nan
                
                
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
                                 metric = drought_metrics_two_threshold(mod_vec=data[:,i,j], lib_path=lib_path,
                                                          perc_onset=perc_onset, perc_termination=perc_termination,
                                                          monthly=monthly, obs_vec=control_ref[:,i,j],
                                                          return_all_tsteps=return_all_tsteps, scale=scale[s],
                                                          add_metrics==(['rel_intensity', 'threshold', 
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
                    data_tim  = ncfile.createVariable('timing',   'i4',('time','lat','lon'), fill_value=miss_val)

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
                    os.system("nccompress -o " + out_file) 


