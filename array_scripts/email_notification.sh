#!/bin/bash

#SBATCH -n 1 #Number of cores 
#SBATCH -t 00:30:00 #Runtime in HH:MM:SS
#SBATCH -p serial_requeue #Partition to submit to 
#SBATCH --mem-per-cpu=64 #Memory per cpu in MB (see also --mem) 
#SBATCH --mail-type=END #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@fas.harvard.edu #Email to which notifications will be sent

exit 0

