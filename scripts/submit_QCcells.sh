#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_cellQC.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_cellQC.err
#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c12


singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif Rscript /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/qc_cells.R --matrix /mnt/scratchb/jmlab/morgan02/Covid/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz --barcodes /mnt/scratchb/jmlab/morgan02/Covid/quant.dir/SIGAA1_cells.txt --sparsitythreshold 0.99 --umithreshold 1000 --mtthreshold 2 --plots /mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA1_QC.pdf --out /mnt/scratchb/jmlab/morgan02/Covid/QC/SIGAA1_barcode-whitelist.txt
