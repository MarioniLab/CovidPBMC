#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_BCR_SIGAB1.out  
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_BCR_SIGAB1.err
#SBATCH --mem=60000
#SBATCH -n1
#SBATCH -c24


singularity exec /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif cellranger vdj --id=SIGAB1_BCR --fastqs=data/fastq --reference=genome_refs/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 --sample=SIGAB1_b --localcores=24 --localmem=60
