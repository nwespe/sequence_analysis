# sequence_analysis
Scripts for running sequence analysis

This is code I wrote for use of John Koschwanez’s mutantanalysis scripts in high throughput on the Odyssey cluster. It consists of three scripts to be run in sequence. All three scripts take as input (1) a directory where input files (fastq, pileup, bam or varscan) are located, (2) an output directory, and (3) a csv file containing names of the samples to be analyzed, with ancestors, clones and pools listed in separate columns. See example csv files for formatting. Each script first checks whether the files it will create already exist and whether the input files it needs are in the directory provided. If the script executes successfully, each script ends with the option to have an email notification sent using the email_notification.sh script. To use this, copy the email_notification.sh script from the common directory to your own directory, open the script and change the email address in line 8 to your own. Take note of the full path to your email script.


The first script finds fastq files for the samples listed in the csv file and creates pileup files. All output files are placed into a new folder with the sample name in a new directory "alignment_files" in the given output directory. The pileup files are created by calling the pairedf2p function from the segtools script. This is a modified version of John’s original segtools script. The pileup functions for all samples are run in parallel through the cluster.   

input: python /full/path/to/script1.py -i /your/path/to/fastq_files -d /your/path/to/output_directory -f /your/path/to/csv_file.csv

output: /your/path/to/output_directory/alignment_files/sample1/sample1.pileup for each individual sample listed in csv_file
The output folder will also contain other files created in pileup process (trimmed.fastq, .intervals, .bam, etc)


The second script finds pileup files for the samples listed in the csv file and creates snp and indel files for each ancestor-clone pair and for each pool (vs. ref). All output files are placed into the relevant subdirectory of a new directory "varscan_files." The varscan files are created by calling the varscan_clone or varscan_pool function from the segtools_modified script. These are submitted as two separate array jobs to the cluster: one for ancestor-clone pairs and one for pools. The input for this script is the same as the first script.

input: python /full/path/to/script1.py -i /your/path/to/pileup_files -d /your/path/to/output_directory -f /your/path/to/csv_file.csv

If -d is omitted, the directory is assumed to be the same one specified by -i.

output: /your/path/to/output_directory/varscan_files/anc_vs_clone/anc_vs_clone.snp, anc_vs_clone.indel
/your/path/to/output_directory/varscan_files/pool/pool.snp, pool.indel


The third script finds bam files (created by script 1) and snp and indel files (created by script 2) and submits these to John’s mutantanalysis script to obtain the comparison output for each ancestor-clone-pool or ancestor-clone(-clone...) group. 
input: python /full/path/to/script3.py -b /your/path/to/bam_files -v /your/path/to/varscan_files -d /your/path/to/output_directory -f /your/path/to/csv_file.csv -m mut_percent -s seg_percent

If -v and/or -d are omitted, the directories are assumed to be the same one specified by -b.

-m and -s are integers between 0-100 that indicate the percent of reads needed to call a mutation (m) or segregant (s). Defaults are m=90 and s=70. Recommend m=35 for diploid clones. 

The output for each comparison is placed into the relevant subdirectory of a new directory "analysis_output-p-m%-s%" where p = pipeline (sawc or cmc), m% = mutation percent, s% = segregation percent (e.g. analysis_output-sawc-m90-s70)
