#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Extract_ADTexprs.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Extract_ADTexprs.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
SCEfiles=$(eval 'echo "$BASE_DIR"/SCE/Covid_SCE-SIGAA1_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA2_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA3_Ab.RDS')
#echo $SCEfiles

DonorFiles=$(eval 'echo "$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv')
#echo $DonorFiles

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/extract_counts.R --SCElist $SCEfiles --DonorList $DonorFiles  --output $BASE_DIR/MOFA/Covid_ADT
