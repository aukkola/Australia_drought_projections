#!/bin/bash

#PBS -P oq98
#PBS -l walltime=02:30:00
#PBS -l mem=5GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/wj02

module load cdo

#Calculate mean annual temperature over land

path="/g/data/w97/amu561/CABLE_AWRA_comparison"

met_path="/g/data/wj02/COMPLIANT_PUBLISHED/HMINPUT/output/AUS-5/BoM/"

#Experiments to process
declare -a experiments=('historical' 'rcp45' 'rcp85')

#output path
outdir=$path'/Global_mean_temp/CABLE_AWRA/'

mkdir $outdir


#Use CABLE/AWRA data to mask out oceans
#Ideally would use models' own land masks but these are not available for all
#models. 
awap_mask=$path/landsea/landsea_mask_AWRA_fixed.nc


#Loop through experiments
for exp in "${experiments[@]}"
do
  
  #List models
  bc_methods=`ls $path/../Steven_CABLE_runs/drought_metrics/3-month/`
  
  for b in $bc_methods
  do
  
    models=`ls $path/../Steven_CABLE_runs/drought_metrics/3-month/${b}/`
  
    for mod in $models
    do
      
      #Outfile
      outdir_mod=$outdir"/"$exp"/"$b"/"$mod
      mkdir -p $outdir_mod
            
      outfile=$outdir_mod"/"${exp}"_"${mod}"_"${b}"_global_mean_temp.nc"

      #If exists, skip model/bc combo
      if [ -e $outfile ];
      then
        echo "$exp, $mod, $b exists, skipping"
        continue
      fi  
        
      #The BC method names are longer in the BoM folder
      #Match using the shorter bc name (CCAM has two folders, need to pick correct one)
      if [ $b == "CCAM" ];
      then
        pattern="*${b}*ISIMIP*"
      else
        pattern="r240x120-${b}-AWAP"
      fi
      
      
      #Find file
      tasmin_file=`ls ${met_path}/${mod}/${exp}/r1i1p1/${pattern}/latest/day/tasmin/*.nc`
      tasmax_file=`ls ${met_path}/${mod}/${exp}/r1i1p1/${pattern}/latest/day/tasmin/*.nc`
      
      #First calculate monthly means from daily data
      tmin_temp="${outdir_mod}/tmin_temp.nc"
      tmax_temp="${outdir_mod}/tmax_temp.nc"

      cdo monmean $tasmin_file $tmin_temp
      cdo monmean $tasmax_file $tmax_temp
      

      #Then calculate mean temperature
      temp_file="${outdir_mod}/temp.nc"
      cdo ensmean $tmin_temp $tmax_temp $temp_file

      #Calculate global mean (excludes Antarctica through GRUN masking)
      cdo -L fldmean $temp_file $outfile
          
      rm $temp_file $tmin_temp $tmax_temp
      
    done #ensembles
    
  done #models
  
done #experiments







