#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Classify_cells_randomforests.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Classify_cells_randomforests.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/classify_cells.R
