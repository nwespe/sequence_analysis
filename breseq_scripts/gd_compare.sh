#!/bin/bash

#SBATCH -n 1 #Number of cores
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem)

script_directory=$1
reference=$2
shift 2
comp_info=("$@")
num_comparisons=${#comp_info[@]}

echo "Received ${num_comparisons} sets of gd files to compare"

tuple=${comp_info[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a comparison <<< "$tuple"  # (title, output_dir, num_samples, gdfile1, [gdfile2,...])
title="${comparison[0]}"
output="${comparison[1]}"
num_samples="${comparison[2]}"
gdfiles_array=("${comparison[3]}")  # start array of file names with first gdfile

i=1  # index of second gdfile in all_files array
while [ $i -lt $num_samples ]  # put all gdfiles into all_files array
do
    comp_index=$[$i+3]
    gdfiles_array+=("${comparison[$comp_index]")
    i=$[$i+1]
done
all_files=$( IFS=' '; echo "${gdfiles_array[*]}" )  # convert array to space-delimited string

cd $output
echo "Now running gdtools COMPARE for $title in $output"
echo "Reference file is $reference"
gdtools COMPARE -o $title.html -r $reference [${all_files}]

exit 0

