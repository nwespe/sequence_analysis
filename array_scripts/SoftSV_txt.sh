#!/bin/bash

#SBATCH -n 1 #Number of cores
#SBATCH -t 3-00:00 #Runtime in D-HH:MM
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH --mem-per-cpu=20480 #Memory per cpu in MB (see also --mem)

samples=("$@")
num_samples=${#samples[@]}

echo "Received bam files and output directories for ${num_samples} samples"

tuple=${samples[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a SAMPLE <<< "$tuple"
BAM="${SAMPLE[0]}"
BASENAME="${BAM##*/}"
NAME="${BASENAME%.bam}"
OUTPUT="${SAMPLE[1]}"
echo "Now running SoftSV for $NAME"
SoftSV -i "${BAM}" -o "${OUTPUT}"

exit 0