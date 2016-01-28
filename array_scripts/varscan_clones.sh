#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=20480 #Memory per cpu in MB (see also --mem) 

script_directory=$1
shift
paired_pileups=("$@")
num_pairs=${#paired_pileups[@]}

echo "Received ${num_pairs} ancestor-clone pairs"

tuple=${paired_pileups[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a PAIR <<< "$tuple"
ANC_PILEUP="${PAIR[0]}"
CLONE_PILEUP="${PAIR[1]}"
NAME="${PAIR[2]}"
OUTPUT="${PAIR[3]}"
cd $OUTPUT
echo "Now running varscan for $NAME in $OUTPUT"
$script_directory/segtools varscan_pair "${ANC_PILEUP}" "${CLONE_PILEUP}" "${NAME}"

exit 0

