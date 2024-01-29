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


#Experiments to process
declare -a experiments=('historical' 'rcp45' 'rcp85' 'rcp26')

#output path
outdir=$path'/Global_mean_temp/CORDEX/'

mkdir $outdir


#Use CABLE/AWRA data to mask out oceans
#Ideally would use models' own land masks but these are not available for all
#models. 
awap_mask=$path/landsea/landsea_mask_AWRA_fixed.nc


#Loop through experiments
for exp in "${experiments[@]}"
do
  
  #List GCMs
  gcm_models=`ls $path"/CORDEX_data/Processed_CORDEX_data/"$exp"/tas/"`
  
  
  for mod in $gcm_models
  do
    
    #List RCMs
    rcm_models=`ls $path"/CORDEX_data/Processed_CORDEX_data/"$exp"/tas/"$mod`

    for rcm in $rcm_models
    do
      
      
      #List ensembles
      ensembles=`ls $path"/CORDEX_data/Processed_CORDEX_data/"$exp"/tas/"$mod"/"$rcm`
      
      for ens in $ensembles
      do
        
        #Find file
        file=`ls ${path}/CORDEX_data/Processed_CORDEX_data/${exp}/tas/${mod}/${rcm}/${ens}/*_Aus.nc`
        
        #Drop path from file name, used to create output file
        file_nopath=`echo "$file" | sed 's!.*/!!'`

        #Outfile
        outdir_mod=$outdir"/"$exp"/"$mod"/"$rcm"/"$ens
        mkdir -p $outdir_mod
        
        outfile=$outdir_mod"/"${file_nopath%"Aus.nc"}"global_mean.nc"


        #Use GRUN file to mask. First resample it to the model resolution
        awap_temp=$outdir_mod/awap_temp.nc
        cdo -L remapnn,$file $awap_mask $awap_temp
        
        
        #Temp file for masking
        temp_mask=$outdir_mod"/temp_mask.nc"
        
        #Mask oceans where max runoff zero
        cdo -L div $file -gtc,0 $awap_temp $temp_mask
        
        #Calculate global mean (excludes Antarctica through GRUN masking)
        cdo -L fldmean $temp_mask $outfile
            
        rm $temp_mask $awap_temp
        
      done #ensembles
    done #RCMs
  done #GCMs
done #experiments







