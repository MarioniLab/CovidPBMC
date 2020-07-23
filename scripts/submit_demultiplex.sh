#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_demux.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_demux.err
#SBATCH --mem=100000
#SBATCH -n1 
#SBATCH -c12

singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif vireo --vartrixData /mnt/scratchb/jmlab/morgan02/Covid/genotypes/vartrix/SIGAA1_alt_matrix.mtx,/mnt/scratchb/jmlab/morgan02/Covid/genotypes/vartrix/SIGAA1_ref_matrix.mtx,/mnt/scratchb/jmlab/morgan02/Covid/quant.dir/SIGAA1_cells.txt --nDonor 4 --outDir /mnt/scratchb/jmlab/morgan02/Covid/demultiplexed/vartrix_vireo/SIGAA1 --randSeed 42 --noPlot
