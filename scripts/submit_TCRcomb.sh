#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/TCR_QC.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/TCR_QC.err
#SBATCH --mem=36000
#SBATCH -n1 
#SBATCH -c1

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/combine_TCR_GEX.R
