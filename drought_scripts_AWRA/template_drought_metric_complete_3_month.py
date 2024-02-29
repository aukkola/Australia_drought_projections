## Adapted code from Anna Ukkola 

#################################
### IMPORT NECESSARY PACKAGES ###
#################################

from netCDF4 import Dataset,num2date 
import numpy as np
import glob
import sys 
import os
import datetime
import xarray as xr

#################
### Set paths ###
#################

lib_path  = "/g/data/w97/amu561/Steven_CABLE_runs/scripts/drought_scripts/functions"

##### ALTERED #####
#data_path = "/g/data/w97/amu561/Steven_CABLE_runs/"
data_path = "/g/data/wj02/COMPLIANT_PUBLISHED/"

scratch_path = '/scratch/w97/amu561/'


# Add lib_path to os directory
sys.path.append(os.path.abspath(lib_path))

#import all droughtmetric calcuation functions
from drought_metrics import *

########################
### DEFINE VARIABLES ###
########################

### Bias Correction ###
bias_corr=str(sys.argv[1])

### Set Models ###
model=str(sys.argv[2])

### Set scenarios ###
scenario=str(sys.argv[3])

### Set variable ###
variable=str(sys.argv[4])

### Set drought metric conditions ###
return_all_tsteps=True

### Set percentile for drought threshold ###
perc=15
##########################################
##########################################
##########################################
##########################################
### Set scale for month aggregation ###
scale=3
##########################################
##########################################
##########################################
##########################################


### Set to monthly threshold calc ###
monthly=True

### Set additional reference data to false ###
obs_ref = False
obs_var = variable  

### Set historical refernce period ###
baseline=[1970,2005]

##########################
### FILE PREPROCESSING ###
##########################

### Define file locations ###
input_variable=['pr']
output_variable=['qtot','sm_total','sm_root']


# Any output data from CCAM has a slightly differenet path
ccam_add=""
if bias_corr=="CCAM":
    ccam_add="CSIRO-CCAM-r3355-"

awra_add=""
if variable in output_variable:
    awra_add='AWRALv6-1-'

if variable in input_variable: # and bias_corr in compliant:
    data_path_var= data_path + "/HMINPUT/output/AUS-5/BoM"

if variable in output_variable: # and bias_corr in compliant:
    data_path_var= data_path + "/HMOUTPUT/output/AUS-5/BoM"



### Create temporary directory ###
temp_dir_path = f"/scratch/w97/amu561/temp/{bias_corr}{model}{scenario}{variable}{str(scale)}"
os.system("mkdir -p " + temp_dir_path)

os.system("mkdir -p " + scratch_path + "/monthly_sums/")


### Get total soil moisture ###

#Need to get s0 (surface, 0-10cm), ss (shallow soil, 10-100cm) and sd (deep soil (100-600cm))
if variable in ["sm_total", "sm_root"]:

    #Variables to fetch
    extra_vars = ["s0", "ss"]
    
    #Add deep soil if calculating total soil moisture
    if variable=="sm_total":
        extra_vars.append("sd")
    
    
    if bias_corr=="CCAM":
        add="ISIMIP2b"
    else:
        add = bias_corr

    for s in range(len(extra_vars)):
        
        print("extra var: ", extra_vars[s])
        #output file name
        sm_out_file = str(scratch_path + "/monthly_sums/" + bias_corr + "_" + model + "_" + 
                          scenario + "_" )
        

        #If output doesn't already exist, create it
        if not os.path.isfile(str(sm_out_file + extra_vars[s] + ".nc")):

            ### Define file locations ###
            sm_files=glob.glob(str(data_path_var + '/' + awra_add + model + '/historical/' + 
                                        '/r1i1p1/' + ccam_add + 'r240x120-ISIMIP2b-AWAP/' + 
                                        'latest/day/' + extra_vars[s] + '/*1960*.nc'))

            sm_files_fut=glob.glob(str(data_path_var + '/' + awra_add + model + '/' + scenario + 
                                              '/r1i1p1/' + ccam_add + 'r240x120-ISIMIP2b-AWAP/' + 
                                              'latest/day/' + extra_vars[s] + '/*2006*.nc'))

            sm_files.extend(sm_files_fut)
        
            print("SM files: ", *sm_files)
        
            temp_sm_file=str(temp_dir_path + "/temp_" + extra_vars[s] + "_file.nc")
            
            files_sm_merge=" ".join(sm_files)

            ### Merge the monthly sum data ###
            os.system("cdo mergetime " + files_sm_merge + " " + temp_sm_file)
            
            #Take monthly mean
            os.system("cdo monmean " + temp_sm_file + " " + sm_out_file + extra_vars[s] + ".nc")

            #Remove temp file
            #os.system("rm " + temp_sm_file)



#Precip, runoff
else:
    
    ### Get historical and future simulations ###

    #historical
    if bias_corr=="CCAM":
        files_hist=str(data_path_var + '/' + awra_add + model + '/historical/r1i1p1/' +
                           ccam_add + 'r240x120-ISIMIP2b-AWAP/latest/day/' +
                           variable + '/*1960*.nc')
    else:
        files_hist=str(data_path_var + '/' + awra_add + model +'/historical/r1i1p1/' +
                           ccam_add + 'r240x120-' + bias_corr + '-AWAP/latest/day/' +
                           variable + '/*1960*.nc')

    #Files to merge                    
    files_all=glob.glob(files_hist)

    #Future
    if bias_corr=="CCAM":
        files_fut=str(data_path_var + '/' + awra_add + model + '/' + scenario + 
                        '/r1i1p1/' + ccam_add + 'r240x120-ISIMIP2b-AWAP/' + 
                        'latest/day/' + variable + '/*2006*.nc')
    else:
        files_fut=str(data_path_var + '/' + awra_add + model + '/' + scenario + 
                           '/r1i1p1/' + ccam_add + 'r240x120-' + bias_corr + 
                           '-AWAP/latest/day/' + variable + '/*2006*.nc')

    files_to_merge_fut=glob.glob(files_fut)

    files_all.extend(files_to_merge_fut)

    ### Location of output file ###
    merged_file = str(scratch_path + "/monthly_sums/" + bias_corr + "_" + 
                      model + "_" + scenario + "_" + variable + ".nc")


    #Some duplicate time steps, need to skip these when merging
    #os.system("export SKIP_SAME_TIME=1")


    #If output file doesn't alreade exist, create it. Else skip 
    if not os.path.isfile(merged_file):

        
        temp_file=str(temp_dir_path + "/temp_" + variable + "_file.nc")
        
        files_merge=" ".join(files_all)
        ### Merge the monthly sum data ###
        os.system("cdo mergetime " + files_merge + " " + temp_file)
        
        #Calculate monthly mean/sum
        #Average for soil moisture, otherwise sum
        #if variable in ["sm_total", "sm_root"] :
        #    os.system("cdo monmean " + temp_file + " " + merged_file)
        #else:
        os.system("cdo monsum " + temp_file + " " + merged_file)
        
        #Remove temp file
        os.system("rm " + temp_file)


        # ds = xr.open_mfdataset(files_all)
        # 
        # #test=ds['s0'].groupby("time.month").sum()
        # test1=ds['s0'].resample(time="1MS").mean(dim="time")
        # 
        # #Rename data variable
        # test1.name = variable
        # 
        # #Write to file
        # test1.to_netcdf(merged_file, format='NETCDF4', 
        #                encoding={variable:{
        #                          'shuffle':True,
        #                          'chunksizes':[12, 681, 40],
        #                          'zlib':True,
        #                          'complevel':5}
        #                          })




        # #Soil moisture
        # if variable == "sm" :
        # 
        #     monthly_s0 = 
        # 
        # 
        #     #Need to sum s0, ss and sd to get total soil moisture 
        #     extra_vars = ["ss", "sd"]
        # 
        #     if bias_corr=="CCAM":
        #         add="ISIMIP2b"
        #     else:
        #         add = bias_corr
        # 
        #     for s in range(len(extra_vars)):
        # 
        #         #file name (leave unfinished to complete below as also need to use it flexibly
        #         #when reading data back in)
        #         sm_out_file = str(scratch_path + "/monthly_sums/" + bias_corr + "_" + model + "_" + 
        #                           scenario + "_")
                





### Delete temporary directory ###
#os.system("rm -r " + temp_dir_path) 

#################
### Load data ###
#################

#Model data
if variable in ["sm_total", "sm_root"]:

    #s0 data (surface)
    fh        = Dataset(str(sm_out_file + "s0.nc"), mode='r')
    s0_data   = fh.variables["s0"][:]

    #ss data (shallow)
    ds_ss     = Dataset(str(sm_out_file + "ss.nc"), mode='r')
    ss_data   = ds_ss.variables["ss"][:]
    
    #Sum s0 and ss to get root zone soil moisture
    all_data  = s0_data.data + ss_data.data 
    
    #Add deep layer if using total soil moisture
    if variable == "sm_total":
        #sd data (deep)
        ds_sd     = Dataset(str(sm_out_file + "sd.nc"), mode='r')
        sd_data   = ds_sd.variables["sd"][:]

        all_data = all_data + sd_data.data
    
    
    data      = all_data
    fh_time   = fh.variables["time"]
    mask      = s0_data.mask


else:
    
    fh       = Dataset(merged_file, mode='r')
    all_data = fh.variables[variable][:] #[yr_ind]
    
    data     = all_data.data 
    mask     = all_data.mask
    fh_time  = fh.variables["time"]

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

miss_val = -999
data[mask==True] = miss_val

control_ref = data

###################################################
### Create output file name and check if exists ###
###################################################
     
### Define output path ###
out_path = f'/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics/{scale}-month/{bias_corr}/{model}'

if not os.path.exists(out_path):    
    os.makedirs(out_path)
##########################################
# ##########################################
# ##########################################
# ##########################################       
#Create output file name
#if variable=='s0':
#    var_name='sm'
#else:
#    var_name=variable
out_file = str(out_path + "/drought_metrics_" + bias_corr + "_" + model + "_" + 
               variable + "_" + scenario + "_" + str(scale) + ".nc")
               
##########################################
##########################################
##########################################
##########################################
##########################################


#############################
### Find baseline indices ###
#############################

#Get dates
ref_years = fh_years
if obs_ref:
    ref_years = obs_years
    

subset = range(np.where(ref_years == baseline[0])[0][0],
                np.where(ref_years == baseline[1])[0][-1] + 1) #Stupid python indexing
                
################################
### Initialise output arrays ###
################################

if return_all_tsteps:
    save_len = len(data)
else:
    save_len = int(len(data)*(perc/100)*2)

duration          = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
rel_intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
rel_intensity_mon = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan

#intensity     = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan
timing        = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    
#tseries       = np.zeros((save_len, len(lat), len(lon))) + miss_val # * np.nan    

if monthly:
    threshold    = np.zeros((12, len(lat), len(lon))) + miss_val # * np.nan
else:
    threshold    = np.zeros((len(lat), len(lon))) + miss_val # * np.nan

#########################
### Calculate metrics ###
#########################

#Loop through grid cells
for i in range(len(lat)):
            
    for j in range(len(lon)):
    
            #Calculate metrics if cell not missing
            if any(~mask[:,i,j]):  
            
                #Calculate metrics
                metric = drought_metrics(mod_vec=data[:,i,j], lib_path=lib_path, perc=perc, 
                                        monthly=monthly, obs_vec=control_ref[:,i,j],
                                        return_all_tsteps=return_all_tsteps, scale=scale,
                                        add_metrics=(['rel_intensity', 'threshold', 'rel_intensity_monthly',
                                        'timing']),
                                        subset=subset, miss_val=miss_val)
        
                ### Write metrics to variables ###
                duration[range(np.size(metric['duration'])),i,j]   = metric['duration']  #total drought duration (months)

                rel_intensity[range(np.size(metric['rel_intensity'])),i,j] = metric['rel_intensity'] #average magnitude

                rel_intensity_mon[range(np.size(metric['rel_intensity_monthly'])),i,j] = metric['rel_intensity_monthly'] #average magnitude
            
                #intensity[range(np.size(metric['intensity'])),i,j] = metric['intensity'] #average intensity
    
                timing[range(np.size(metric['timing'])),i,j]       = metric['timing']    #drought timing (month index)

                #tseries[range(np.size(metric['tseries'])),i,j]       = metric['tseries']    #drought timing (month index)

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
time      = ncfile.createVariable("time", 'i4', ('time',))

if monthly:
    month = ncfile.createVariable("month", 'i4', ('month',))

#Create data variables
data_dur  = ncfile.createVariable('duration', 'f8',('time','lat','lon'), fill_value=miss_val)
data_mag  = ncfile.createVariable('rel_intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_rel  = ncfile.createVariable('rel_intensity_by_month','f8',('time','lat','lon'), fill_value=miss_val)

#data_int  = ncfile.createVariable('intensity','f8',('time','lat','lon'), fill_value=miss_val)
data_tim  = ncfile.createVariable('timing',   'i4',('time','lat','lon'), fill_value=miss_val)
#data_ts   = ncfile.createVariable('tseries',   'i4',('time','lat','lon'), fill_value=miss_val)

#Create data variable for threshold
if monthly:
    data_thr = ncfile.createVariable('threshold', 'f8',('month','lat','lon'), fill_value=miss_val)
else:
    data_thr = ncfile.createVariable('threshold', 'f8',('lat','lon'), fill_value=miss_val)


#Set variable attributes
longitude.units = 'degrees_east'
latitude.units  = 'degrees_north'
time.units      = 'days since 1960-01-01'

time.calendar   = 'gregorian'

data_dur.long_name = 'drought event duration (no. months)'
data_mag.long_name = 'drought event relative intensity (%)'
data_rel.long_name = 'drought month relative intensity (%)'

#data_int.long_name = 'drought event intensity (mm)'
data_tim.long_name = 'drought event timing (binary drought/non-drought index)'
data_thr.long_name = 'drought threshold (mm)'
#data_ts.long_name  = 'original time series'

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

#data_int[:,:,:] = intensity
data_tim[:,:,:] = timing
#data_ts[:,:,:]  = tseries

if monthly:    
    data_thr[:,:,:] = threshold
else:
    data_thr[:,:] = threshold

# Close the file
ncfile.close()


#Compress file
os.system("nccompress " + out_file)
os.system("mv " + out_path + "/tmp.nc_compress/*" + variable + "*.nc " + out_path)
#os.system("rm -r " + out_path + "/tmp.nc_compress")
# jobs=`qselect -u amu561`
# for J in $jobs
# do
# echo $J
# if [[ "$J" =~ "109307765.gadi-pbs" ]]; then
# echo "skip"
# else
# qdel $J
# fi
# done













