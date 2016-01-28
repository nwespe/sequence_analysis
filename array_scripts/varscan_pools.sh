#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem) 

script_directory=$1
shift
pool_pileups=("$@")
num_pools=${#pool_pileups[@]}

echo "Received ${num_pools} pools"

tuple=${pool_pileups[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a POOL <<< "$tuple"
PILEUP="${POOL[0]}"
BASENAME="${PILEUP##*/}"
NAME="${BASENAME%.pileup}"
OUTPUT="${POOL[1]}"
echo "Now running varscan for $NAME"
$script_directory/segtools varscan_pool "${PILEUP}" "${NAME}" "${OUTPUT}"

exit 0

