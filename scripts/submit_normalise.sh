#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_normalise.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_normalise.err
#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c1

singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif Rscript /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/norm_data.R --matrixlist /mnt/scratchb/jmlab/morgan02/Covid/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz --barcodeslist /mnt/scratchb/jmlab/morgan02/Covid/quant.dir/SIGAA1_cells.txt         --whitelists /mnt/scratchb/jmlab/morgan02/Covid/QC/SIGAA1_barcode-whitelist.txt --sparsity 0.95 --output /mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS --sizefactors /mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SumFactors.tsv --plots /mnt/scratchb/jmlab/morgan02/Covid/reports/Covid_norm.pdf --force --outputqc
