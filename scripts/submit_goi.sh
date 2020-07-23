#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_goi.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_goi.err
#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c1

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/fast_goi.R
