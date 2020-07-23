#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Train_MOFA_singlegroup.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Train_MOFA_singlegroup.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1

## Run MOFA treating all cells as from the same sample - does it pick out any batch effects?

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
SCEfiles=$(eval 'echo "$BASE_DIR"/MOFA/Covid_ADT_longCPM.txt.gz')
#echo $SCEfiles

DonorFiles=$(eval 'echo "$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv')
#echo $DonorFiles

LogFile=$(echo "$BASE_DIR"/logs/MOFA_singleGroup_TrainPy.log)

# turn off HDF5 file locking
export HDF5_USE_FILE_LOCKING='FALSE'

python /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/train_mofa.py --counts-list=$SCEfiles --donor-list=$DonorFiles --factors=50 --runmode=medium --format=long --ignore-groups --MOFAoutput=$BASE_DIR/MOFA/Covid_ADT_singleGroup-MOFA.hdf5 --output-prefix=$BASE_DIR/MOFA/Covid_ADT_singleGroup --log=$LogFile
