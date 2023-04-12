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

met_path="/g/data/ia39/australian-climate-service/release/CORDEX-CMIP6/output/AUS-15/BOM/"

#Experiments to process
declare -a experiments=('historical' 'ssp126' 'ssp370')

#output path
outdir=$path'/Global_mean_temp/BARPA/'

mkdir $outdir


#Use monthly processed mrro data to mask out oceans
#Note this is stored on scratch so will cease to exist after some time
#can be recreated using the drought metric codes
mrro_mask="/scratch/w97/amu561/monthly_sums_BARPA/historical_EC-Earth-Consortium-EC-Earth3_r1i1p1f1_mrro.nc"


#Loop through experiments
for exp in "${experiments[@]}"
do
  
  #list models
  models=`ls $met_path`

  
    for mod in $models
    do
      
      #Skip ECMWF-ERA5, only historical evaluation data
      if [ $mod == "ECMWF-ERA5" ];
      then
        continue
      fi
      
      
      #List ensemble members
      ensembles=`ls ${met_path}/${mod}/${exp}/`

      #Loop through ensemble members
      for ens in $ensembles
      do
      
        
        #Outfile
        outdir_mod=$outdir"/"$exp"/"$mod"/"${ens}
        mkdir -p $outdir_mod
              
        outfile=$outdir_mod"/"${exp}"_"${mod}"_"${ens}"_global_mean_temp.nc"

        #If exists, skip model/bc combo
        if [ -e $outfile ];
        then
          echo "$exp, $mod, $b exists, skipping"
          continue
        fi  
        
        
        #Find files
        #Some mean temp files seem to be in a folder called "tas" and others
        #in "tasmean". Need to deal with this
        
        if [ -d "${met_path}/${mod}/${exp}/${ens}/BOM-BARPA-R/v1/mon/tas/" ];
        then
          tas_files=`ls ${met_path}/${mod}/${exp}/${ens}/BOM-BARPA-R/v1/mon/tas/*.nc`
        else
          tas_files=`ls ${met_path}/${mod}/${exp}/${ens}/BOM-BARPA-R/v1/mon/tasmean/*.nc`
        fi      
        

        #Then calculate mean temperature
        temp_file=${outdir_mod}"/temp.nc"
        cdo -L mergetime $tas_files $temp_file

        #Calculate global mean (excludes Antarctica through GRUN masking)
        cdo -L fldmean $temp_file $outfile
          
        rm $temp_file
      
    done #ensembles
    
  done #models
  
done #experiments







