#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=6000 #Memory per cpu in MB (see also --mem) 

while getopts :d:r: opt
do
    case $opt in
        d) script_directory=$OPTARG
            ;;
        r) ref_seq=$OPTARG
            ;;
    esac
done
shift $((OPTIND-1))

samples=("$@")

tuple=${samples[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a sample <<< "$tuple"
input="${sample[0]}"
output="${sample[1]}"

#echo "${sample[@]}"
echo "Now running picard-tools CollectAlignmentSummaryMetrics for $input."

java -jar $PICARD/CollectAlignmentSummaryMetrics.jar R=$ref_seq I=$input O=$output VALIDATION_STRINGENCY=LENIENT

exit 0
