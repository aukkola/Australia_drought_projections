#!/bin/bash

#PBS -P oq98
#PBS -q normalbw
#PBS -l storage=gdata/er4+scratch/w97+gdata/hh5+gdata/wj02+gdata/eg3+gdata/w97+gdata/oq98
#PBS -l walltime=03:00:00
#PBS -l ncpus=1
#PBS -l mem=60gb
#PBS -j oe


module use /g/data/hh5/public/modules
module load conda/analysis3

python3 /g/data/oq98/amu561/Steven_CABLE_runs/scripts/drought_scripts_CABLE/template_drought_metric_complete_3_month.py \
$scenario $model $variable $co2 $out_file $scale

#"MODEL" "VARIABLE" "CO2" "OUTFILE"
