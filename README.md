# sequence_analysis
Scripts for running sequence analysis

This is code I wrote for use of John Koschwanez’s mutantanalysis scripts in high throughput on the Odyssey cluster. It consists of three scripts to be run in sequence. All three scripts take as input (1) a directory where your fastq files are located and (2) a csv file containing names of the samples to be analyzed, with ancestors, clones and pools listed in separate columns. See example csv files for formatting. Each script first checks whether the files it will create already exist and whether the input files it needs are in the directory provided. If the script executes successfully, each script ends with the option to have an email notification sent using the email_notification.sh script. To use this, copy the email_notification.sh script from the common directory to your own directory, open the script and change the email address in line 8 to your own. Take note of the full path to your email script.

The first script finds fastq files for the samples listed in the csv file, sorts paired fastq files into individual subdirectories under a new directory "fastq_pairs," and creates pileup files. All output files are placed into the relevant subdirectory of fastq_pairs. The pileup files are created by calling the pairedf2p function from the segtools_modified script. This is John’s original segtools script with separate varscan functions added. The pileup functions for all samples are run in parallel through the cluster.   
input: python /full/path/to/script1.py /your/path/to/fastq_files /your/path/to/csv_file.csv
output: /your/path/to/fastq_files/fastq_pairs/sample1/sample1.pileup for each individual sample listed in csv_file
The output folder will also contain other files created in pileup process (.sam, .bam, etc)

The second script finds pileup files for the samples listed in the csv file and creates snp and indel files for each ancestor-clone pair and for each pool (vs. ref). All output files are placed into the relevant subdirectory of a new directory "varscan_files." The varscan files are created by calling the varscan_clone or varscan_pool function from the segtools_modified script. These are submitted as two separate array jobs to the cluster: one for ancestor-clone pairs and one for pools. The input for this script is the same as the first script.
input: python /full/path/to/script1.py /your/path/to/fastq_files /your/path/to/csv_file.csv
output: /your/path/to/fastq_files/varscan_files/anc_vs_clone/anc_vs_clone.snp, anc_vs_clone.indel
/your/path/to/fastq_files/varscan_files/pool/pool.snp, pool.indel

The third script finds bam files (created by script 1) and snp and indel files (created by script 2) and submits these to John’s mutantanalysis script to obtain the comparison output for each ancestor-clone-pool or ancestor-clone(-clone...) group. 

The input for this script has more parameters than the other scripts and uses "flags" to indicate them:
input: python /full/path/to/script3.py -d /your/path/to/fastq_files -f /your/path/to/csv_file.csv -p pipeline -m mut_percent -s seg_percent

Required:
-d is the directory and -f is the csv file used for scripts 1 and 2

Optional:
-p is the pipeline you wish to run: either 'sawc' (segregant analysis with clone) or 'cmc' (compare multiple clones). The pipeline run by the program is determined by the format of csv file, so this parameter is optional. For sawc, label columns Ancestor, Clone, Pool_1, Pool_2, etc. - see example file 1 for appropriate formatting to create the necessary files. For cmc, label columns Clone_1, Clone_2, etc. - see example file 2 for formatting. The appropriate varscan files will be created if the same csv file is used for script 2.
-m and -s are integers between 0-100 that indicate the percent of reads needed to call a mutation (m) or segregant (s). Defaults are m=90 and s=70. Recommend m=35 or 40 for diploid clones. 

The output for each comparison is placed into the relevant subdirectory of a new directory "analysis_output-p-m%-s%" where p = pipeline chosen, m% = mutation percent, s% = segregation percent (e.g. analysis_output-sawc-m70-s70)

A fourth script (copynumber.py) runs the varscan copynumber function for ancestor-clone pairs. The output files (anc_vs_clone.copynumber) is placed into the relevant varscan_files/ subfolder, the same place as the .snp and .indel files. I am currently working on a method to visualize these data.
