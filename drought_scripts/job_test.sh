#!/bin/bash

#PBS -P oq98
#PBS -q normal
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/eg3+gdata/w97+gdata/oq98
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=10gb
#PBS -j oe


module use /g/data/hh5/public/modules
module load conda/analysis3

mkdir -p /scratch/w97/amu561/test

echo "done here"

ls /scratch/w97/amu561/test

python3 /g/data/w97/amu561/Steven_CABLE_runs/drought_scripts/template_drought_metric_complete_3_month.py "CCAM" "CNRM-CERFACS-CNRM-CM5" "rcp45" "pr"
