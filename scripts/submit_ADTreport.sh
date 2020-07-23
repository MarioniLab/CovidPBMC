#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_ADTreport.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_ADTreport.err
#SBATCH --mem=24000
#SBATCH -n1 
#SBATCH -c1


Rscript /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/ADT_report.R --CITEdirectory /mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAA1_CITE --ADTSCE /mnt/scratchb/jmlab/morgan02/Covid/SCE/Covid_SCE.RDS --plots /mnt/scratchb/jmlab/morgan02/Covid/reports/ADT-QC.pdf
