#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_TCR_SIGAC1.out  
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_TCR_SIGAC1.err
#SBATCH --mem=60000
#SBATCH -n1
#SBATCH -c24


singularity exec /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif cellranger vdj --id=SIGAC1_TCR --fastqs=data/fastq --reference=genome_refs/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 --sample=SIGAC1_t --localcores=24 --localmem=60
