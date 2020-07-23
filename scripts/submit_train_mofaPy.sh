#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Train_MOFA.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Train_MOFA.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
SCEfiles=$(eval 'echo "$BASE_DIR"/MOFA/Covid_ADT_SIGAA1_cpm.txt.gz,"$BASE_DIR"/MOFA/Covid_ADT_SIGAA2_cpm.txt.gz,"$BASE_DIR"/MOFA/Covid_ADT_SIGAA3_cpm.txt.gz')
#echo $SCEfiles

DonorFiles=$(eval 'echo "$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv')
#echo $DonorFiles

LogFile=$(echo "$BASE_DIR"/logs/MOFA_TrainPy.log)

# turn off HDF5 file locking
export HDF5_USE_FILE_LOCKING='FALSE'

python /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/train_mofa.py --counts-list=$SCEfiles --donor-list=$DonorFiles --factors=50 --format=wide  --runmode=medium --MOFAoutput=$BASE_DIR/MOFA/Covid_ADT_MOFA.hdf5 --output-prefix=$BASE_DIR/MOFA/Covid_ADT --log=$LogFile
