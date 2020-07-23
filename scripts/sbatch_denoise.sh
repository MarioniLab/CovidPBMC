#! /usr/bin/bash

#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c1

# 1) UMI threshold
# 2) Output plots
# 3) Output prefix

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
SCEfiles=$(eval 'echo "$BASE_DIR"/SCE/Covid_SCE-SIGAA1_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA2_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA3_Ab.RDS')

DonorFiles=$(eval 'echo "$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv')

matrixFiles=$(eval 'echo "$BASE_DIR"/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA2_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA3_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz')

cellFiles=$(eval 'echo "$BASE_DIR"/quant.dir/SIGAA1_cells.txt,"$BASE_DIR"/quant.dir/SIGAA2_cells.txt,"$BASE_DIR"/quant.dir/SIGAA3_cells.txt')


JOB="Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/denoise_ADT.R --SCElist $SCEfiles --DonorList $DonorFiles --matrixlist $matrixFiles --whitelist $cellFiles --IgControls IgG1_Ctrl,IgG2a_Ctrl,IgG2b_Ctrl,IgG2b_RatCtrl --plots $2 --output $3 --UMIthreshold $1"

eval $JOB

