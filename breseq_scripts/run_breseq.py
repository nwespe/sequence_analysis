#!/usr/bin/env python
# -*- coding: utf-8 -*-
# run_breseq.py

""" Description: match sample names to fastq files and run breseq.
    User input will be a directory where fastq files are present
    and a csv file with sample names.
"""

__author__= "Nichole Wespe"

import subprocess
import os
import sys
import csv
import fnmatch
import argparse

script_directory = os.path.abspath(os.path.dirname(sys.argv[0]))

parser = argparse.ArgumentParser(description='Tell me where to find fastq files, where to put the output files, '
                                             'and what samples to analyze.')
parser.add_argument('-i', '--input', required=True, help='Full path of directory containing fastq files. Files can be '
                                                         'in subdirectories of this directory.')
parser.add_argument('-o', '--output', required=True, help='Full path of directory for file output. Must be different '
                                                          'from input directory containing fastq files.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing sample names. Names must '
                                                            'be in fastq file names.')
parser.add_argument('-r', '--reference', required=False, default=os.path.join(script_directory, 'reference.fa'),
                    help='Full path to reference file. Default is script_directory/reference.fa')

args = parser.parse_args()

fastq_directory = args.input
directory = args.output
csv_file = args.csv_file
ref_file = args.reference


def input_eval(fastq_directory, directory):
    
    """ Determine if inputs are a directory and a readable csv file
    """
    try:
        os.path.isdir(fastq_directory)
    except:
        print directory + " is not a directory"
        sys.exit(1)

    try:
        os.path.isdir(directory)
    except:
        print directory + " is not a directory"
        sys.exit(1)

    #if fastq_directory == directory:
     #   print 'Will not output to same directory as fastq files. Please specify a different output directory.'
     #   sys.exit(1)
        
    print "Now working on files in " + directory


def read_csv_file(csv_file):
    
    """ Read csv file of samples to analyze and create list of sample names.
    """
    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)

    # make set of sample names to pass to other functions
    individual_samples = set()  # set of samples to analyze
    for row in csv_reader:
        for sample in row:
            if sample: individual_samples.add(sample)

    print "Received directions to analyze " + str(len(individual_samples)) + " unique samples"

    return individual_samples


def check_output(samples):

    """ Check for breseq output files in output directory.
    """
    gd_files = []
    for r, d, f in os.walk(directory):
        for item in f:
            if fnmatch.fnmatch(item, '*output.gd'):
                gd_files.append(os.path.abspath(os.path.join(r, item)))

    gd_not_found = set()
    gd_found = set()

    for sample in samples:
        name_pattern = '*'+sample+'*'
        matched_file = fnmatch.filter(gd_files, name_pattern)  # generates list of strings
        if len(matched_file) == 0: gd_not_found.add(sample)  # not found
        else: gd_found.add(sample)  # found and unique
        
    if bool(gd_found):
        print 'One or more GenomeDiff files found in output folder for: ' + ', '.join(str(f) for f in gd_found) + \
              '. Now running breseq analysis for remaining samples.'

    return gd_not_found


def find_fastq_files(samples): # input is set of entries from csv file for which no pileups were found

    """ Matches pairs of fastq files, moves pair into new subdirectory if necessary,
        and returns list of tuples (R1 file, R2 file, sample name) with full pathnames.
    """
    # make list of all fastq files in given directory tree
    fastq_files = []
    for r, d, f in os.walk(fastq_directory):
        for item in f:
            if item[0] != '.' and fnmatch.fnmatch(item, '*.fastq*' or '*.fq*'):
                fastq_files.append(os.path.abspath(os.path.join(r, item)))
  
    # find fastq files corresponding to samples (input list)
    sample_file_dict = {}
    for sample in samples:
        sample_file_dict[sample] = fnmatch.filter(fastq_files, '*'+sample+'*')
        # creates dict containing the sample name and a sublist with one or two (or more, which is bad) files

    fastq_pairs = []
    fastq_single = []
    not_unique = []
    not_found = []
    unpaired = []
    for key, file_list in sample_file_dict.iteritems():  # dict is sample name = fastq list
        if len(file_list) == 0: not_found.append(key)
        elif len(file_list) >2: not_unique.append(key)
        elif len(file_list) == 1:  # check if unpaired or single read
            fastq_file = file_list[0] 
            if fnmatch.fnmatch(fastq_file, '*R[12].f*'): unpaired.append(fastq_file)
            else:
                output_dir = os.path.abspath(os.path.join(directory, 'breseq_analysis', key))
                if not os.path.exists(output_dir): os.makedirs(output_dir)
                fastq_single.append((fastq_file, key, output_dir))
        elif len(file_list) == 2: 
            R1 = fnmatch.filter(file_list, '*R1.f*')
            R2 = fnmatch.filter(file_list, '*R2.f*')
            if len(R1) == 1 and len(R2) == 1:  # if both files are found and are unique
                R1_file = R1[0]
                R2_file = R2[0]
                output_dir = os.path.abspath(os.path.join(directory, 'breseq_analysis', key))
                if not os.path.exists(output_dir): os.makedirs(output_dir)
                fastq_pairs.append((R1_file, R2_file, key, output_dir))
            else: not_unique.append(key)  # if the 2 files in list are not an R1-R2 pair

    if bool(not_unique): print 'ID names not unique for: ' + ', '.join(str(f) for f in not_unique)
    if bool(not_found): print 'fastq files not found for: ' + ', '.join(str(f) for f in not_found)
    if bool(unpaired): print 'Only one fastq file found for paired read: ' + ', '.join(str(f) for f in unpaired) 

    return fastq_single, fastq_pairs


def fastq_to_breseq(fastq_info):

    """ fastq_info is a list of tuples, either (file, sample_name, output_dir) or
        (R1, R2, sample_name, output_dir)
    """
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    array_num = len(fastq_info)-1
    bash_file = os.path.join(script_directory, 'fastq_to_breseq.sh')  # 'array_scripts',
    output = os.path.join(diag_dir, 'breseq_%A-%a.out')
    error = os.path.join(diag_dir, 'breseq_%A-%a.err')
    arg_list = []
    if len(fastq_info[0]) == 3:  # fastq_info is for single reads
        single_or_paired = 'single'
        print "Now running breseq analysis for " + str(len(fastq_info)) + " fastq files"
        for single in fastq_info:
            arg_list.append(' '.join(single))
    else:
        single_or_paired = 'paired'
        print "Now running breseq analysis for " + str(len(fastq_info)) + " pairs of fastq files"
        for pair in fastq_info:
            arg_list.append(' '.join(pair))
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output),
                    bash_file, script_directory, ref_file, single_or_paired]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch fastq_to_breseq.sh executed, return code is " + str(returncode)
    return returncode


def main():  # run analyses of fastq files using functions defined above

    """ check input,
        extract data from csv file
        check if GenomeDiff files (part of breseq output) exist,
        run breseq if gd_not_found - in array format
    """
    # initialize lists
    fastq_single = []
    fastq_pairs = []
    gd_not_found = set()

    # set default exit statuses
    run_status_single = 0; run_status_pairs = 0
    
    input_eval(fastq_directory, directory)  # check input
    individual_samples = read_csv_file(csv_file)

    # empty sets or lists initialized above evaluate to False; if any have members created by previous function,
    # "if" statement below is True and dependent function is executed
    gd_not_found = check_output(individual_samples)

    if not gd_not_found:
        print "Found breseq GenomeDiff files in output directory for all samples."
        sys.exit(0)
    else: fastq_single, fastq_pairs = find_fastq_files(gd_not_found)

    #fastq_single, fastq_pairs = find_fastq_files(individual_samples)

    if fastq_single:
        run_status_single = fastq_to_breseq(fastq_single)

    if fastq_pairs: 
        run_status_pairs = fastq_to_breseq(fastq_pairs)

    if run_status_single or run_status_pairs != 0:
        # error while executing breseq function; status is 0 if function not performed
        print "Error encountered while executing breseq function. Check sbatch output files."
        sys.exit(1)
    else:
        print ("To receive an email when analyses have been completed, confirm that your email address is in line 8 of "
               "a copy of script email_notification.sh. Then type the following into the command line, substituting in "
               "the sbatch job ID and full path to script:")
        print "sbatch --dependency=afterok:<jobid> </your/path/to/email_notification.sh>"
        print "To check job status, type squeue -u <username> -j<jobid> into the command line."
        print "When finished, use run_gdtools.py script to compare multiple breseq analyses if desired."
        sys.exit(0)


if __name__ == "__main__": 
    main()
