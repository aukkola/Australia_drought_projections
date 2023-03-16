#!/bin/bash

#PBS -P oq98
#PBS -l walltime=00:30:00
#PBS -l mem=5GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l storage=gdata/w97

module load cdo

#Calculate mean annual temperature over land

path="/g/data/w97/amu561/CABLE_AWRA_comparison"

#Paths to AGCD data
tmax_path="/g/data/zv2/agcd/v1/tmax/mean/r005/01month/"
tmin_path="/g/data/zv2/agcd/v1/tmin/mean/r005/01month/"

#output directory
outdir=$path/Observed_tair/Obs_mean_time_series
mkdir -p $outdir


#Merge tmax and min files

tmax_files=`find $tmax_path/*.nc`
cdo mergetime $tmax_files $outdir/tmax_merged.nc

tmin_files=`find $tmin_path/*.nc`
cdo mergetime $tmin_files $outdir/tmin_merged.nc


#Calculate mean temperature
cdo ensmean $outdir/tmax_merged.nc $outdir/tmin_merged.nc $outdir/tmean.nc


#Calculate spatial time series
cdo fldmean $outdir/tmean.nc $outdir/Monthly_mean_temperature_time_series_Aus.nc

#Remove temp files
rm $outdir/*_merged.nc $outdir/tmean.nc


