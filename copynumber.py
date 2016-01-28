#!/usr/bin/env python
# -*- coding: utf-8 -*-
# copynumber.py

""" Description: Create copynumber files using VarScan.
    User input will be a directory where fastq files are present
    and a csv file describing analyses to be done. All comparisons
    to ancestor should be specified as clones, not pools.
"""

__author__= "Nichole Collins"

import subprocess
import os
import sys
import csv
import re
import fnmatch
import argparse

script_directory = os.path.dirname(sys.argv[0])

parser = argparse.ArgumentParser(description = "Specify directory and csv file.")
parser.add_argument('-d', '--directory', required=True, help='Full path of directory containing fastq files. Files can be in subdirectories of this directory.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing desired comparisons and names of ancestors, clones, pools. See sample documents for formatting.')
args = parser.parse_args()

directory = args.directory
csv_file = args.csv_file


def input_eval(directory, csv_file):

    """ Determine if inputs are a directory and a readable csv file
    """
    try:
        os.path.isdir(directory)
    except:
        print directory + " is not a directory"
        sys.exit(1)
##    try:
##        csv_open = open(csv_file, 'rU')
##        csv_reader = csv.reader(csv_open)
##        #some other function here to check file?
##    except:
##        print csv_file + " is not readable csv file."
##        return 1

    print "Now working on files in " + directory


def read_csv_file(csv_file):
    """ Read csv file of analyses to be done and create input for find_pileups
        function. First row must be header row with "Ancestor, Clone_x"
        Sends one set of ancestor-clone pairs.
    """
# problem: can't read csv file with uneven numbers of entries in rows
# hack: format all empty cells as text in excel before saving as csv

    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)
    row1 = next(csv_reader)

    # get number of pools
    clone_regex = re.compile('(Clone).*')
    clones = ([c.group(0) for cell in row1 for c in [clone_regex.search(cell)] if c])
    num_clones = len(clones)
    pool_regex = re.compile('(Pool).*')
    pools = ([p.group(0) for cell in row1 for p in [pool_regex.search(cell)] if p])
    if pools: print "Use csv file specifying all comparisons to ancestor as clones."

    # get indices for ancestor, clone, pool(s)
    ancestor_index = row1.index('Ancestor')
    clone_indices = []
    for c in clones:
        clone_indices.append(row1.index(c))

    # make sets to pass to other functions
    individual_samples = set() # set of samples to make pileups of
    ancestor_clone_pairs = set() # pairs for varscan analysis
    for row in csv_reader:
        ancestor = row[ancestor_index]
        individual_samples.add(ancestor)
        list_of_clones = []
        for c in clone_indices:
            clone = row[c]
            if clone:
                individual_samples.add(clone)
                ancestor_clone_pairs.add((ancestor, clone))
                list_of_clones.append(clone)
 
    print "Received directions to analyze " + str(len(individual_samples)) + \
          " unique samples and " + str(len(ancestor_clone_pairs)) + " anc-clone pairs. " \

    return individual_samples, ancestor_clone_pairs


def check_cn_files(pairs):

    """ Check for snp and indel files for ancestor_clone_pairs
    """
    # find all snp and indel files in directory tree:
    cn_files = []
    for r, d, f in os.walk(directory):
        for item in f:
            if fnmatch.fnmatch(item, '*.copynumber'):
                cn_files.append(os.path.abspath(os.path.join(r, item)))

    if not cn_files:  # no files yet
        samples = set() # initialize set of individual samples to find pileups for
        for pair in pairs:
            for sample in pair: samples.add(sample)
        print 'Now finding pileup files for ' + str(len(samples)) + ' samples'
        return samples, pairs

    # below script should only run if any copynumber files were found

    # make dictionary of (key=anc_vs_clone, value=(ancestor, clone))
    # anc_vs_clone name is used to search copynumber files
    anc_clone_dict = {}
    for pair in pairs: # ancestor_clone_pairs is a set of tuples
        anc_clone_name = str(pair[0] + '_vs_' + pair[1]) # name of analysis
        anc_clone_dict[anc_clone_name] = pair

    # create sets of pairs, pools not already analyzed
    filtered_pairs = set()
    analyzed_pairs = set()
    for key, pair in anc_clone_dict.iteritems():
        key_pattern = '*'+key+'*'
        matched_cn = fnmatch.filter(cn_files, key_pattern)
        if not matched_cn: filtered_pairs.add(pair) # missing file(s)
        else: analyzed_pairs.add(pair) # both files exist

    filtered_samples = set() # initialize set of individual samples to find pileups for
    for pair in filtered_pairs:
        for sample in pair: filtered_samples.add(sample)

    return filtered_samples, filtered_pairs


def find_pileups(samples, pairs):

    """ Get set of ancestor-clone pairs and find pileup files in filepath.
    """

    pileup_files = []
    for r, d, f in os.walk(directory):
        for item in f:
            if fnmatch.fnmatch(item, '*.pileup'):
                pileup_files.append(os.path.abspath(os.path.join(r, item)))

    pileups_not_found = set()
    pileups_not_unique = set()

    pileups_dict = {}
    for sample in samples:
        name_pattern = '*'+sample+'*'
        matched_file = fnmatch.filter(pileup_files, name_pattern) #generates list of strings
        if len(matched_file) == 1: pileups_dict[sample] = matched_file[0]
        elif len(matched_file) == 0: pileups_not_found.add(sample) #not found
        elif len(matched_file) >1: pileups_not_unique.add(sample) #not unique

    paired_pileups = []
    for pair in pairs:
        pileup_pair = [] #list of ancestor and clone
        for item in pair:
            if item in pileups_dict: item_pileup = pileups_dict[item]; pileup_pair.append(item_pileup) #adds both ancestor and clone to one list
        if len(pileup_pair) == 2:
            anc_clone_name = str(pair[0] + '_vs_' + pair[1]) #makes name for comparison
            output_dir = os.path.abspath(os.path.join(directory, 'varscan_files', anc_clone_name))
            if not os.path.exists(output_dir): os.makedirs(output_dir)
            paired_pileups.append((pileup_pair[0], pileup_pair[1], anc_clone_name, output_dir))

    if pileups_not_found: print 'Pileups not found for: ' + ', '.join(str(f) for f in pileups_not_found)
    if pileups_not_unique: print 'ID names not unique for: ' + ', '.join(str(f) for f in pileups_not_unique)

    return paired_pileups, pileups_not_found #, pileups_not_unique


def batch_varscan_copynumber(paired_pileups): # paired_pileups is a list of tuples (anc.pileup, clone.pileup, anc_clone_name, output_dir)

    """ Run varscan copynumber routine for ancestor-clone pairs from find_pileups
        function.
    """
    array_num = len(paired_pileups)-1
    bash_file = os.path.join(script_directory, 'array_scripts', 'copynumber.sh')
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    output = os.path.join(diag_dir, 'varscan_copynumber_%A-%a.out')
    error = os.path.join(diag_dir, 'varscan_copynumber_%A-%a.err')
    print "Now creating varscan files for " + str(len(paired_pileups)) + " ancestor-clone pairs"
    arg_list = []
    for pair in paired_pileups:
        arg_list.append(' '.join(pair))
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file, script_directory]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch copynumber.sh executed, return code is " + str(returncode)
    return returncode


def main(): # run analyses of pileup files using functions defined above

    """ check input,
        extract data from csv file
        check if copynumber files exist
        run varscan function if copynumer files don't exist - in array format
    """

    # initialize lists
    samples = set(); pairs = set()
    fastq_pairs = []
    pileups_not_found = set()

    # set default exit statuses
    var_cn_status = 0; new_cn_file = 0

    input_eval(directory, csv_file) # check input
    individual_samples, ancestor_clone_pairs = read_csv_file(csv_file)

    # empty sets or lists initialized above evaluate to False; if any have members created by previous function, "if" statement below is True and dependent function is executed
    samples, pairs = check_cn_files(ancestor_clone_pairs)

    if not samples:
        print "All varscan copynumber analyses have been performed."
        sys.exit(0)

    paired_pileups, pileups_not_found = find_pileups(samples, pairs)

    if pileups_not_found:
        print "Not all pileup files were found. Run script 1 to generate pileups."
        sys.exit(1) # pileup files need to be created

    # now all required pileup files should exist - create snp, indel files if needed
    if paired_pileups:
        new_cn_file = 1
        var_cn_status = batch_varscan_copynumber(paired_pileups)
    if var_cn_status !=0: # error while making varscan files; status is 0 if function not performed
        print "Error encountered while making copynumber files. Check error log for varscan jobs."
        sys.exit(1)

    if new_cn_file == 1: # if varscan function was run
        print "Varscan analyses being run."
        print "To receive an email when varscan files have been completed, confirm that your email address is in line 8 of a copy of script email_notification.sh. Then type the following into the command line, substituting in the sbatch job ID(s) and full path to script:"
        print "sbatch --dependency=afterok:<jobid1>:<jobid2> </your/path/to/email_notification.sh>"
        print "To check job status, type squeue -u <username> -j<jobid> into the command line."
        sys.exit(0)
    else:
        print "Error encountered while executing varscan function."
        sys.exit(1)


if __name__ == "__main__":
    main()




