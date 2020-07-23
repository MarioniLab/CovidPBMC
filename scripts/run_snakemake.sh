#!/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Snakemake_run.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Snakemake_run.err
#SBATCH --mem=12000
#SBATCH -n1
#SBATCH -c1

DT=$(date '+%d-%m_%H-%M')
BASEDIR=/mnt/scratchb/jmlab/morgan02/Covid/
ERR=$BASEDIR/logs/Snakemake_$DT.e
OUT=$BASEDIR/logs/Snakemake_$DT.o

cd $BASEDIR

snakemake -s CovidPBMC/snakemake_files/Snakefile_CRUK --cluster "sbatch -c {cluster.nCPUs} --nodes={cluster.nodes} --job-name={cluster.name} --mem={cluster.memory} --output={cluster.output} --error={cluster.error} " --cluster-config $BASEDIR/CovidPBMC/cluster_config/CRUK_cluster.json --jobs=100 --latency-wait=250 --verbose --rerun-incomplete -rp  --until normalise_gex # -R
