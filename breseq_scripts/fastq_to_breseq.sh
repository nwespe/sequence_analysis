#!/bin/bash

#SBATCH -n 1 #Number of cores
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem)

script_directory=$1
reference=$2
single_or_paired=$3
shift 3
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
    echo "Now running breseq for $NAME in $OUTPUT"
    echo "Reference file is $reference"
    breseq -p -n "$NAME" -r $reference "$READ"
elif [ $single_or_paired = 'paired' ]
then
    READ1="${SAMPLE[0]}"
    READ2="${SAMPLE[1]}"
    NAME="${SAMPLE[2]}"
    OUTPUT="${SAMPLE[3]}"
    cd $OUTPUT
    echo "Now running pairedf2p for $NAME in $OUTPUT"
    echo "Reference file is $reference"
    breseq -p -n "$NAME" -r $reference "$READ1" "$READ2"
fi

exit 0
