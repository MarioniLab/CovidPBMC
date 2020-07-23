#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_vartrix.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_vartrix.err
#SBATCH --mem=148000
#SBATCH -n1 
#SBATCH -c12

singularity exec  -B /mnt/scratchb  /mnt/scratchb/jmlab/morgan02/Covid/singularity/CovidPBMC.sif vartrix --vcf /mnt/scratchb/jmlab/morgan02/Covid/genotypes/genomeProject/hg38/genome1K.phase3.SNP_AF5e2.chr1toX.withoutChr.vcf.gz --bam /mnt/scratchb/jmlab/morgan02/Covid/SIGAA1_CITE/outs/possorted_genome_bam.bam --fasta /mnt/scratchb/jmlab/morgan02/Covid/genome_refs/hg38_SarsCov2/fasta/genome.fa --cell-barcodes /mnt/scratchb/jmlab/morgan02/Covid/quant.dir/SIGAA1_cells.txt --threads 20 --log-level info --scoring-method coverage --out-matrix /mnt/scratchb/jmlab/morgan02/Covid/genotypes/vartrix/SIGAA1_alt_matrix.mtx --ref-matrix /mnt/scratchb/jmlab/morgan02/Covid/genotypes/vartrix/SIGAA1_ref_matrix.mtx
