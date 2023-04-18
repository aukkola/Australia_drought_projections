#!/bin/bash

#PBS -P oq98
#PBS -q normal
#PBS -j oe
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/w97+gdata/oq98
#PBS -l walltime=00:10:00
#PBS -l ncpus=1
#PBS -l mem=1gb

path="/g/data/oq98/amu561/CABLE_AWRA_comparison/"

PBS_PATH=${path}/scripts/drought_scripts_CABLE/pbs_jobs

mkdir -p ${PBS_PATH}

cd ${path}/scripts/drought_scripts_CABLE

### Set scale ###
scale=3

#Loop through models
for model in CSIRO-BOM-ACCESS1-0 CNRM-CERFACS-CNRM-CM5 NOAA-GFDL-GFDL-ESM2M MIROC-MIROC5; do

  #Loop through variables
  for variable in qtot sm; do
    
    #Loop through CO2 options
    for co2 in CO2 noCO2; do

      #Loop through scenarios
      for scenario in rcp45 rcp85; do


        #No data for constant CO2 for RCP8.5 scenario, skip
        if [ $co2 == "noCO2" -a $scenario == "rcp85" ];
        then
          continue
        fi
        
        
        #Output directory
        outdir="/g/data/w97/amu561/Steven_CABLE_runs/drought_metrics_CABLE/${scale}-month/${co2}/${model}/"  
        mkdir -p $outdir
        
        #Output file
        out_file=${outdir}/drought_metrics_CABLE_${co2}_${model}_${variable}_${scenario}.nc 
             
        #If output doesn't exist, process
        if [ ! -f ${out_file} ]; then

          echo "here"
          pbs_job=${PBS_PATH}/job_${model}_${co2}_${variable}_${scenario}.sh
          
          echo $scenario $model $variable $co2 $out_file $pbs_job
          
          cp job_template_${scale}.sh $pbs_job
          # wait 
          # sed -i "s/MODEL/${model}/g" $pbs_job
          # wait
          # sed -i "s/VARIABLE/${variable}/g" $pbs_job
          # wait
          # sed -i "s/CO2/${co2}/g" $pbs_job
          # wait
          # sed -i "s/OUTFILE/${out_file}/g" $pbs_job
          # wait
          qsub -v "scenario=$scenario","model=$model","variable=$variable","co2=$co2","out_file=$out_file","scale=$scale" $pbs_job
          
        fi 
      done 
    done  
  done
done
