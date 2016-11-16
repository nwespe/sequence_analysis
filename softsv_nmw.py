#!/usr/bin/env python
# -*- coding: utf-8 -*-
# softsv_nmw.py

""" Description
    This script will run SoftSV for bam files
    User input will be a directory where bam files are present
    and a csv file describing samples to be analyzed.
"""

__author__= "Nichole Wespe"

import subprocess
import os
import sys
import csv
import re
import fnmatch
import argparse

script_directory = os.path.dirname(sys.argv[0])

parser = argparse.ArgumentParser(description = "Specify a directory where bam files are present and a csv file "
                                               "describing samples to be analyzed.")
parser.add_argument('-b', '--bam_directory', required=True, help='Full path of directory containing bam files. Files can '
                                                             'be in subdirectories of this directory.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing names of samples.')
parser.add_argument('-d', '--directory', default='same', help='Full path of directory for file output. If omitted, '
                                                              'output directory is same as csv file directory.')
args = parser.parse_args()

bam_directory = args.bam_directory
csv_file = args.csv_file
directory = args.directory
if directory is 'same':
    directory = os.path.abspath(os.path.dirname(csv_file))

def input_eval():

    """ Determine if input is a directory
    """
    try:
        os.path.isdir(bam_directory)
    except:
        print bam_directory + " is not a directory"
        sys.exit(1)

    print "Now working on files in " + bam_directory


def read_csv_file():  # TODO: this currently finds any entry in a csv file past row 1

    """ Read csv file of sample names
    """

    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)
    row1 = next(csv_reader)

    # make set of sample names to pass to other functions
    individual_samples = set() # set of samples to analyze bam files
    for row in csv_reader:
        for sample in row:
            if sample: individual_samples.add(sample)

    print "Received directions to analyze " + str(len(individual_samples)) + " unique samples"

    return individual_samples


def filter_softsv_files(samples):

    # TODO: find pre-existing report files, need to know what output is
    data_files = []
    for r, d, f in os.walk(directory):  # checks provided output directory
        for item in f:
            if fnmatch.fnmatch(item, '*softsv.txt'):  ## TODO: common part of the name of the SoftSV output file
                data_files.append(os.path.abspath(os.path.join(r, item)))

    # match data files to samples
    filtered_samples = []
    analyzed_samples = []
    for sample in samples:
        sample_file = fnmatch.filter(data_files, '*'+sample+'*')  ## TODO: assumes sample is part of output file name
        if sample_file: analyzed_samples.append(sample)
        else: filtered_samples.append(sample)

    if analyzed_samples:
        print "Found existing SoftSV output file for " + ', '.join(analyzed_samples)

    return filtered_samples


def find_bam_files(samples):

    # find bam files
    bam_files = []
    for r, d, f in os.walk(bam_directory):
        for item in f:
            if fnmatch.fnmatch(item, '*final.bam'):  ## TODO: uses the final bam file created by the current pipeline
                bam_files.append(os.path.abspath(os.path.join(r, item)))
    print "Found " + str(len(bam_files)) + "final bam files"

    # match bam files to samples
    matched_samples = []
    non_unique_samples = []
    for sample in samples:
        sample_bam = fnmatch.filter(bam_files, '*'+sample+'*')
        if len(sample_bam) > 1: non_unique_samples.append(sample)
        else: matched_samples.append((sample, sample_bam[0]))

    if non_unique_samples:
        print "More than one final bam file found for " + ', '.join(non_unique_samples)
    print "Matched " + str(len(matched_samples)) + " samples to bam files"
    return matched_samples


def run_SoftSV(matched_samples):

    # TODO: modify for arguments for SoftSV command
    # ref_seq = os.path.join(script_directory, 'mutantanalysis', 'ref_seq', 's288c_sgd.fa')
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    output = os.path.join(diag_dir, 'softsv_report_%A-%a.out')
    error = os.path.join(diag_dir, 'softsv_report_%A-%a.err')

    analysis_dir = os.path.join(directory, 'softsv_reports-' + os.path.splitext(os.path.basename(csv_file))[0])
    arg_list = []
    array_num = len(matched_samples)-1
    for sample, bam in matched_samples:
        output_dir = os.path.abspath(os.path.join(analysis_dir, sample))  # put reports in same dir as csv
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        #output_file = os.path.abspath(os.path.join(analysis_dir, sample+'_softsv.txt'))
        sample_args = (bam, output_dir)  # create tuple containing args for SoftSV command
        arg_list.append(' '.join(sample_args))  # convert tuple to space-delimited string and add to list

    bash_file = os.path.join(script_directory, 'SoftSV_txt.sh')
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch SoftSV_txt.sh executed, return code is " + str(returncode)
    return returncode


def main():  # run analyses of bam files using functions defined above

    """ check input,
        extract data from csv file
        find bam files
        run SoftSV_txt.sh script
    """
    # initialize lists
    individual_samples = set()

    input_eval() # check input

    individual_samples = read_csv_file()

    filtered_samples = filter_softsv_files(individual_samples)

    matched_samples = find_bam_files(filtered_samples)

    run_SoftSV(matched_samples)

    print "SoftSV is running on " + str(len(matched_samples)) + " samples. To check status, type squeue " \
                                                                "-u <username> -j<jobid> into the command line."
    print ("To receive an email when analyses have been completed, confirm that your email address is in line 8 of a "
           "copy of script email_notification.sh. Then type the following into the command line, substituting in the "
           "sbatch job ID and full path to script:")
    print "sbatch --dependency=afterok:<jobid1> </your/path/to/email_notification.sh>"
    sys.exit(0)


if __name__ == "__main__":
    main()

