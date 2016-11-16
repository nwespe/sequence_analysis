#!/usr/bin/env python
# -*- coding: utf-8 -*-
# run_gdtools.py

""" Description: match sample names to output.gd files and run gdtools COMPARE.
    User input will be a directory where output.gd files are present
    and a csv file with sample names.
"""

__author__= "Nichole Wespe"

import subprocess
import os
import shutil
import sys
import csv
import fnmatch
import argparse

script_directory = os.path.abspath(os.path.dirname(sys.argv[0]))

parser = argparse.ArgumentParser(description='Tell me where to find GenomeDiff files, where to put the output files, '
                                             'and what samples to analyze.')
parser.add_argument('-i', '--input', required=True, help='Full path of directory containing GenomeDiff files. '
                                                         'Files can be in subdirectories of this directory.')
parser.add_argument('-o', '--output', required=True, help='Full path of directory for file output. Must be different '
                                                          'from input directory containing fastq files.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing sample names. Names must '
                                                            'be in fastq file names.')
parser.add_argument('-r', '--reference', required=False, default=os.path.join(script_directory, 'reference.fa'),
                    help='Full path to reference file. Default is script_directory/reference.fa')

args = parser.parse_args()

input_directory = args.input
output_directory = args.output
csv_file = args.csv_file
ref_file = args.reference


def input_eval(input_directory, output_directory):

    """ Determine if inputs are a directory and a readable csv file
    """
    try:
        os.path.isdir(input_directory)
    except:
        print input_directory + " is not a directory"
        sys.exit(1)

    try:
        os.path.isdir(output_directory)
    except:
        print output_directory + " is not a directory"
        sys.exit(1)

    print "Now working on files from " + input_directory


def read_csv_file(csv_file):

    """ Read csv file of samples to analyze and create list of sample names.
    """
    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)

    # make set of comparisons to pass to other functions
    comparisons = set()  # set of comparisons to make
    all_samples = set()
    for row in csv_reader:
        compare = []
        for sample in row:
            if sample:
                all_samples.add(sample)
                compare.append(sample)
        comparisons.add(tuple(compare))  # makes a set of tuples

    print "Received directions to create " + str(len(comparisons)) + " comparisons"

    return all_samples, comparisons


def find_gd_files(all_samples):

    """ Check for breseq output files in output directory.
    """
    gd_files = []
    for r, d, f in os.walk(input_directory):
        for item in f:
            if fnmatch.fnmatch(item, '*output.gd'):
                gd_files.append(os.path.abspath(os.path.join(r, item)))

    gd_not_found = set()
    gd_not_unique = set()
    gd_file_dict = {}

    for sample in all_samples:
        name_pattern = '*'+sample+'*'
        matched_file = fnmatch.filter(gd_files, name_pattern)  # generates list of strings
        if len(matched_file) == 0: gd_not_found.add(sample)  # not found
        elif len(matched_file) > 1: gd_not_unique.add(sample)
        else:
            gd_file_dict[sample] = matched_file[0]

    if bool(gd_not_unique):
        print 'More than one GenomeDiff file found for: ' + ', '.join(str(f) for f in gd_not_unique) + \
              '. Sample excluded from analyses.'
    if bool(gd_not_found):
        print 'GenomeDiff files not found for: ' + ', '.join(str(f) for f in gd_not_found) + \
              '. Now running gd comparisons using remaining samples.'

    return gd_file_dict


def copy_gd_files(comparisons, gd_file_dict):

    """ copy and rename output.gd files for each comparison to new directory, compile information for each comparison
    """
    all_comp_info = []
    for sample_set in comparisons:
        title = '-'.join(sample_set)  # create name for comparison from all sample names
        output_dir = os.path.abspath(os.path.join(output_directory, title))  # create output folder from name
        if not os.path.exists(output_dir): os.makedirs(output_dir)
        comp_info = [title, output_dir, str(len(sample_set))]  # create list of info for comparison
        for sample in sample_set:  # copy all output.gd files to output folder and rename with sample name
            gd_file = gd_file_dict[sample]
            gd_file_copy = os.path.abspath(os.path.join(output_dir, sample+'.gd'))
            shutil.copy2(gd_file, gd_file_copy)
            comp_info.append(gd_file_copy)  # add sample.gd file to list of comparison info
        all_comp_info.append(comp_info)

    return all_comp_info


def gd_compare(all_comp_info):

    """ all_comp_info is a list of lists: [title, output_dir, num_samples, gdfile1, [gdfile2,...]]
    """
    diag_dir = os.path.join(output_directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    array_num = len(all_comp_info)-1
    bash_file = os.path.join(script_directory, 'gd_compare.sh')  # 'array_scripts',
    output = os.path.join(diag_dir, 'gd_compare_%A-%a.out')
    error = os.path.join(diag_dir, 'gd_compare_%A-%a.err')

    print "Now running gdtools compare analysis for " + str(len(all_comp_info)) + " sets of gd files"
    arg_list = []
    for comparison in all_comp_info:
        arg_list.append(' '.join(comparison))  # convert list to space-delimited string and add to list
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output),
                    bash_file, script_directory, ref_file]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch gd_compare.sh executed, return code is " + str(returncode)
    return returncode


def main():  # run analyses of fastq files using functions defined above

    """ check input,
        extract data from csv file
        run gdtools compare in SLURM array format
    """
    # initialize lists

    # set default exit status
    run_status = 0

    input_eval(input_directory, output_directory)  # check input
    individual_samples, comparisons = read_csv_file(csv_file)

    # empty sets or lists initialized above evaluate to False; if any have members created by previous function,
    # "if" statement below is True and dependent function is executed

    gd_file_dict = find_gd_files(individual_samples)

    all_comp_info = copy_gd_files(comparisons, gd_file_dict)

    run_status = gd_compare(all_comp_info)

    if run_status != 0:
        # error while executing breseq function; status is 0 if function not performed
        print "Error encountered while executing gdtools function. Check sbatch output files."
        sys.exit(1)
    else:
        print ("To receive an email when analyses have been completed, confirm that your email address is in line 8 of "
               "a copy of script email_notification.sh. Then type the following into the command line, substituting in "
               "the sbatch job ID and full path to script:")
        print "sbatch --dependency=afterok:<jobid> </your/path/to/email_notification.sh>"
        print "To check job status, type squeue -u <username> -j<jobid> into the command line."
        sys.exit(0)


if __name__ == "__main__":
    main()
