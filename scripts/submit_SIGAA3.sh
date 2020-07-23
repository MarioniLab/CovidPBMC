#! /usr/bin/bash

#SBATCH -e logs/run_CITE_cellranger.sample=SIGAA3.err
#SBATCH -o logs/run_CITE_cellranger.sample=SIGAA3.out
#SBATCH --mem=60000
#SBATCH -n1
#SBATCH -c24
#SBATCH --job-name=run_CITE_cellranger.sample=SIGAA3

cd /mnt/scratchb/jmlab/morgan02/Covid//cellranger_output/; rm -rf SIGAA3_CITE; singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif cellranger count --id=SIGAA3_CITE --transcriptome=/mnt/scratchb/jmlab/morgan02/Covid/genome_refs/hg38_SarsCov2 --libraries=/mnt/scratchb/jmlab/morgan02/Covid/data/lib_files/SIGAA3_5gex_protein_library.csv --feature-ref=/mnt/scratchb/jmlab/morgan02/Covid/data/ref_files/SIGAA3_5gex_protein_feature_ref.csv  --localcores=24 --localmem=60 --nosecondary; cd /mnt/scratchb/jmlab/morgan02/Covid/
