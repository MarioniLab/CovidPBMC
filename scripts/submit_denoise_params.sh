#! /usr/bin/bash

#SBATCH --output=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_denoise_params.out
#SBATCH --error=/mnt/scratchb/jmlab/morgan02/Covid/logs/Run_denoise_params.err
#SBATCH --mem=1000
#SBATCH -n1 
#SBATCH -c1

# Vary the UMI threshold for DSB normalisation.
UMI_t=(5 10 20 50 100)
BASE_DIR="/mnt/scratchb/jmlab/morgan02/Covid"

for umit in ${UMI_t[@]};
do
    echo $umit
    OUTPUT=$(echo "$BASE_DIR/Projections/Covid_ADT_UMI-$umit")
    OUTPLOT=$(echo "$BASE_DIR/reports/Covid_ADT_UMI-$umit")
    
    JOB="sbatch /mnt/scratchb/jmlab/morgan02/Covid/scripts/sbatch_denoise.sh $umit $OUTPLOT $OUTPUT"
    eval $JOB   

done



