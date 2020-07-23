#! /usr/bin/bash

#SBATCH --output=logs/BCR_TCR_report.out
#SBATCH --error=logs/BCR_TCR_report.err
#SBATCH --mem=48000
#SBATCH -n1
#SBATCH -c1

R CMD /mnt/scratchb/jmlab/morgan02/Covid/CovidPBMC/src/TCR_BCR_report.R --logs=/mnt/scratchb/jmlab/morgan02/Covid/logs/TCR_BCR_report.log --plot=/mnt/scratchb/jmlab/morgan02/Covid/reports/TCR_BCR_report.pdf --TCRdirectory=/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAC1_TCR,/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAC2_TCR,/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAC3_TCR --BCRdirectory=/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAB1_BCR,/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAB2_BCR,/mnt/scratchb/jmlab/morgan02/Covid/cellranger_output/SIGAB3_BCR
