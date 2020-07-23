#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/BgShift_all.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/BgShift_all.err
#SBATCH --mem=24000
#SBATCH --time=08:00:00
#SBATCH -n1 
#SBATCH -c1

BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"
MATRIXfiles=$(eval 'echo "$BASE_DIR"/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA2_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz,"$BASE_DIR"/cellranger_output/SIGAA3_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz')
echo $MATRIXfiles

FeatureFiles=$(eval 'echo "$BASE_DIR"/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/features.tsv.gz,"$BASE_DIR"/cellranger_output/SIGAA2_CITE/outs/raw_feature_bc_matrix/features.tsv.gz,"$BASE_DIR"/cellranger_output/SIGAA3_CITE/outs/raw_feature_bc_matrix/features.tsv.gz')
echo $FeatureFiles

BarcodeFiles=$(eval 'echo "$BASE_DIR"/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/barcodes.tsv.gz,"$BASE_DIR"/cellranger_output/SIGAA2_CITE/outs/raw_feature_bc_matrix/barcodes.tsv.gz,"$BASE_DIR"/cellranger_output/SIGAA3_CITE/outs/raw_feature_bc_matrix/barcodes.tsv.gz')
echo $BarcodeFiles

WhiteList=$(eval 'echo "$BASE_DIR"/QC/SIGAA1_barcode-whitelist.txt,"$BASE_DIR"/QC/SIGAA2_barcode-whitelist.txt,"$BASE_DIR"/QC/SIGAA3_barcode-whitelist.txt')
echo $WhiteList

Rscript /mnt/scratchb/jmlab/morgan02/Covid/scripts/bgshift_cpm.R --matrixlist $MATRIXfiles --featurelist $FeatureFiles  --barcodeslist $BarcodeFiles --whitelist $WhiteList  --output $BASE_DIR/MOFA --plots $BASE_DIR/reports/
