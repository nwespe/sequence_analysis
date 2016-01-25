#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem) 

script_directory=$1
single_or_paired=$2
shift 2
fastq_info=("$@")
num_samples=${#fastq_info[@]}

echo "Received ${num_samples} fastq samples"

tuple=${fastq_info[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a SAMPLE <<< "$tuple"
if [ $single_or_paired = 'single' ]
then
    READ="${SAMPLE[0]}"
    NAME="${SAMPLE[1]}"
    OUTPUT="${SAMPLE[2]}"
    cd $OUTPUT
    echo "Now running singlef2p for $NAME in $OUTPUT"
    $script_directory/segtools singlef2p -a "$READ" -c "$NAME" -r $script_directory/mutantanalysis/ref_seq/s288c_sgd.fa
elif [ $single_or_paired = 'paired' ]
then
    READ1="${SAMPLE[0]}"
    READ2="${SAMPLE[1]}"
    NAME="${SAMPLE[2]}"
    OUTPUT="${SAMPLE[3]}"
    cd $OUTPUT
    echo "Now running pairedf2p for $NAME in $OUTPUT"
    $script_directory/segtools pairedf2p -a "$READ1" -b "$READ2" -c "$NAME" -r $script_directory/mutantanalysis/ref_seq/s288c_sgd.fa
fi

exit 0
