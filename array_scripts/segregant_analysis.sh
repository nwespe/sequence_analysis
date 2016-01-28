#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 10:00:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=10240 #Memory per cpu in MB (see also --mem) 

while getopts :c:d:f:m:p:r:s: opt
do
    case $opt in
        c) ref_chrom=$OPTARG
            ;;
        d) script_directory=$OPTARG
            ;;
        f) ref_feat=$OPTARG
            ;;
        m) mut_percent=$OPTARG
            ;;
        p) pipeline=$OPTARG
            ;;
        r) ref_seq=$OPTARG
            ;;
        s) seg_percent=$OPTARG
            ;;
    esac
done
shift $((OPTIND-1))

samples=("$@")


tuple=${samples[${SLURM_ARRAY_TASK_ID}]}
IFS=' ' read -a sample <<< "$tuple"
title="${sample[0]}"
anc_name="${sample[1]}"
anc_bam="${sample[2]}"
output_dir="${sample[3]}"

#echo "${sample[@]}"
echo "Now running $pipeline mutantanalysis for $title."

if [ $pipeline = "sawc" ]
then
    clone_name="${sample[4]}"
    clone_bam="${sample[5]}"
    clone_snp="${sample[6]}"
    clone_indel="${sample[7]}"
    pool_name="${sample[8]}"
    pool_bam="${sample[9]}"
    pool_snp="${sample[10]}"
    pool_indel="${sample[11]}"
    python $script_directory/mutantanalysis/mutantanalysis.py \
           --title $title --output_dir $output_dir --ref_seq $ref_seq \
           --ref_feat $ref_feat --ref_chrom_num $ref_chrom \
           --ancestor_bam $anc_bam --ancestor_name $anc_name \
           --clone_bam $clone_bam --clone_snp $clone_snp \
           --clone_indel $clone_indel --clone_name $clone_name \
           --pool_bam $pool_bam --pool_name $pool_name \
           --pool_snp $pool_snp --pool_indel $pool_indel \
           --mut_percent $mut_percent --seg_percent $seg_percent
elif [ $pipeline = "cmc" ]
then
    num_clones="${sample[4]}"
    declare -a clones_array
    i="1"
    until [ $i -gt $num_clones ]
    do
        name_index=$[$i*4+1]
        bam_index=$[$i*4+2]
        snp_index=$[$i*4+3]
        indel_index=$[$i*4+4]
        clone_info=("('${sample[$name_index]}', '${sample[$bam_index]}', '${sample[$snp_index]}', '${sample[$indel_index]}', $mut_percent)")
        clones_array=("${clones_array[@]}" "${clone_info}")
        i=$[$i+1]
    done
    all_clones=$( IFS=','; echo "${clones_array[*]}" )
    cd $script_directory/mutantanalysis
    python -c "import mutantanalysis; mutantanalysis.compare_multiple_clones('$title', '$anc_name', '$anc_bam', '$ref_seq', '$ref_chrom', '$ref_feat', '$output_dir', [${all_clones}], False)"
    echo "mutantanalysis.compare_multiple_clones('$title', '$anc_name', '$anc_bam', '$ref_seq', '$ref_chrom', '$ref_feat', '$output_dir', [${all_clones}], False)"
fi

exit 0
