#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Blish_randomforest.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Blish_randomforest.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/train_classifiers.R
