#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_CITE_SIGAA1.out  
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/run_CITE_SIGAA1.err
#SBATCH --mem=60000
#SBATCH -n1
#SBATCH -c24

singularity exec /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif cellranger count --id=SIGAA1_CITE --transcriptome=genome_refs/hg38_SarsCov2 --libraries=data/lib_files/SIGAA1_5gex_protein_library.csv --feature-ref=data/ref_files/SIGAA1_5gex_protein_feature_ref.csv  --localcores=24 --localmem=60 --nosecondary
