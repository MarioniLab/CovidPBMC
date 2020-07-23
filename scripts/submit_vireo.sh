#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_emptydrops.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_emptydrops.err
#SBATCH --mem=48000
#SBATCH -n1 
#SBATCH -c12


singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif Rscript /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/call_cells.R --input /mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAA1_CITE/outs/raw_feature_bc_matrix/matrix.mtx.gz --out /mnt/scratchb/jmlab/morgan02/Covid/quant.dir/SIGAA1_cells.txt --fdrthreshold 0.001 --umithreshold 100 --logs /mnt/scratchb/jmlab/morgan02/Covid/reports/SIGAA1_emptyDrops.pdf
