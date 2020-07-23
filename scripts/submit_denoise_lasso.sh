#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_LASSO_denoise.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_LASSO_denoise.err
#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c1

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
SCEfiles=$(eval 'echo "$BASE_DIR"/SCE/Covid_SCE-SIGAA1_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA2_Ab.RDS,"$BASE_DIR"/SCE/Covid_SCE-SIGAA3_Ab.RDS')
#echo $SCEfiles

DonorFiles=$(eval 'echo "$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA1/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA2/donor_ids.tsv,"$BASE_DIR"/demultiplexed/cellSNP_vireo/SIGAA3/donor_ids.tsv')
#echo $DonorFiles

matrixFiles=$(eval 'echo "$BASE_DIR"/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA2_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA3_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz')
#echo $matrixFiles

cellFiles=$(eval 'echo "$BASE_DIR"/quant.dir/SIGAA1_cells.txt,"$BASE_DIR"/quant.dir/SIGAA2_cells.txt,"$BASE_DIR"/quant.dir/SIGAA3_cells.txt')
#echo $cellFiles

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/lasso_denoise.R --SCElist $SCEfiles --DonorList $DonorFiles  --IgControls IgG1_Ctrl,IgG2a_Ctrl,IgG2b_Ctrl,IgG2b_RatCtrl --plots $BASE_DIR/reports/Covid_ADT- --output $BASE_DIR/Projections/Covid_ADT
