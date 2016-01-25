#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sequence_analysis.py

""" Description: organize fastq files and create pileups.
    User input will be a directory where fastq files are present
    and a csv file describing analyses to be done.
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

parser = argparse.ArgumentParser(description = "Tell me where to find fastq files, where to put the output files, and what samples to align.")
parser.add_argument('-i', '--fastq', required=True, help='Full path of directory containing fastq files. Files can be in subdirectories of this directory.')
parser.add_argument('-d', '--directory', required=True, help='Full path of directory for file output.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing names of ancestors, clones, and/or pools. See example documents for formatting.')
args = parser.parse_args()

fastq_directory = args.fastq
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
        function. First row must be header row with "Ancestor, Clone, Pool_x"
        Sends two sets: one of ancestor-clone pair and one of pools.
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
    num_pools = len(pools)

    # get indices for ancestor, clone, pool(s)
    ancestor_index = row1.index('Ancestor')
    clone_indices = []
    for c in clones:
        clone_indices.append(row1.index(c))
    pool_indices = []
    for p in pools:
        pool_indices.append(row1.index(p))

    # make sets to pass to other functions
    individual_samples = set() # set of samples to make pileups of
    ancestor_clone_pairs = set() # pairs for varscan analysis
    anc_mult_clone_groups = []
    set_of_pools = set() # pools for varscan analysis
    acp_trios = set() # trios for final analysis
    if len(clone_indices) >= 1 and not pool_indices: # multiple clone columns and no pools
        pipeline = "cmc"
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
            anc_mult_clone_groups.append((ancestor, list_of_clones))
    elif len(clone_indices) == 1 and pool_indices: # one clone column and pools
        pipeline = "sawc"
        clone_index = clone_indices[0]
        for row in csv_reader:
            ancestor = row[ancestor_index]
            clone = row[clone_index]
            individual_samples.add(ancestor)
            individual_samples.add(clone)
            ancestor_clone_pairs.add((ancestor, clone))
            for p in pool_indices:
                pool = row[p]
                if pool:
                    individual_samples.add(pool)
                    set_of_pools.add(pool)
                    acp_trios.add((ancestor, clone, pool))

    print "Received directions to analyze " + str(len(individual_samples)) + \
          " unique samples, " + str(len(ancestor_clone_pairs)) + " anc-clone pairs, " \
          + str(len(anc_mult_clone_groups)) + " multi-clone groups, " \
          + str(len(set_of_pools)) + " pools, and " + str(len(acp_trios)) + " trios."

    return individual_samples, ancestor_clone_pairs, anc_mult_clone_groups, set_of_pools, acp_trios, pipeline


def check_pileups(samples):

    """ Check for pileup files in filepath.
    """
    
    pileup_files = []
    for r, d, f in os.walk(directory):
        for item in f:
            if fnmatch.fnmatch(item, '*.pileup'):
                pileup_files.append(os.path.abspath(os.path.join(r, item)))

    pileups_not_found = set()
    pileups_not_unique = set()
    pileups_found = set()

    for sample in samples:
        name_pattern = '*'+sample+'*'
        matched_file = fnmatch.filter(pileup_files, name_pattern) #generates list of strings
        if len(matched_file) == 0: pileups_not_found.add(sample) #not found
        elif len(matched_file) == 1: pileups_found.add(sample) #found and unique
        elif len(matched_file) >1: pileups_not_unique.add(sample) #not unique
        
    if bool(pileups_found): print 'Pileup found for: ' + ', '.join(str(f) for f in pileups_found)
    if bool(pileups_not_unique): print 'More than one pileup found for: ' + ', '.join(str(f) for f in pileups_not_unique)

    return pileups_not_found


def find_fastq_files(samples): # input is set of entries from csv file for which no pileups were found

    """ Matches pairs of fastq files, moves pair into new subdirectory if necessary,
        and returns list of tuples (R1 file, R2 file, sample name) with full pathnames.
        TODO: Also return list of single fastq files.
    """
    # make list of all fastq files in given directory tree
    fastq_files = []
    for r, d, f in os.walk(fastq_directory):
        for item in f:
            if item[0] != '.' and fnmatch.fnmatch(item, '*.fastq' or '*.fq'):
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
    for key, file_list in sample_file_dict.iteritems(): # dict is sample name = fastq list
        if len(file_list) == 0: not_found.append(key)
        elif len(file_list) >2: not_unique.append(key)
        elif len(file_list) == 1: # check if unpaired or single read
            fastq_file = file_list[0] 
            if fnmatch.fnmatch(fastq_file, '*R[12].f*'): unpaired.append(fastq_file)
            else:
                output_dir = os.path.abspath(os.path.join(directory, 'alignment_files', key))
                if not os.path.exists(output_dir): os.makedirs(output_dir)
                fastq_single.append((fastq_file, key, output_dir))
        elif len(file_list) == 2: 
            R1 = fnmatch.filter(file_list, '*R1.f*')
            R2 = fnmatch.filter(file_list, '*R2.f*')
            if len(R1) == 1 and len(R2) == 1: #if both files are found and are unique
                R1_file = R1[0]
                R2_file = R2[0]
                output_dir = os.path.abspath(os.path.join(directory, 'alignment_files', key))
                if not os.path.exists(output_dir): os.makedirs(output_dir)
                fastq_pairs.append((R1_file, R2_file, key, output_dir))
            else: not_unique.append(key) # if the 2 files in list are not an R1-R2 pair

    if bool(not_unique): print 'ID names not unique for: ' + ', '.join(str(f) for f in not_unique)
    if bool(not_found): print 'fastq files not found for: ' + ', '.join(str(f) for f in not_found)
    if bool(unpaired): print 'Only one fastq file found for paired read: ' + ', '.join(str(f) for f in unpaired) 
    return fastq_single, fastq_pairs


def fastq_to_pileup(fastq_info): # fastq_info is a list of tuples, either (file, sample_name, output_dir) or (R1, R2, sample_name, output_dir)
    
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    array_num = len(fastq_info)-1
    bash_file = os.path.join(script_directory, 'array_scripts', 'fastq_to_pileup.sh')
    output = os.path.join(diag_dir, 'pileup_%A-%a.out')
    error = os.path.join(diag_dir, 'pileup_%A-%a.err')
    arg_list = []
    if len(fastq_info[0]) == 3: # fastq_info is for single reads
        single_or_paired = 'single'
        print "Now creating pileups for " + str(len(fastq_info)) + " fastq files"
        for single in fastq_info:
            arg_list.append(' '.join(single))
    else:
        single_or_paired = 'paired'
        print "Now creating pileups for " +str(len(fastq_info)) + " pairs of fastq files"
        for pair in fastq_info:
            arg_list.append(' '.join(pair))
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file, script_directory, single_or_paired]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch fastq_to_pileup.sh executed, return code is " + str(returncode)
    return returncode


def main(): # run analyses of fastq files using functions defined above

    """ check input,
        extract data from csv file
        check if pileup files exist,
        make pileups if pileups_not_found - in array format
            - optional: send new pileups created to varscan function?
    """
   
    # initialize lists
    fastq_single = []
    fastq_pairs = []
    pileups_not_found = set()

    # set default exit statuses
    pileup_status_single = 0; pileup_status_pairs = 0
    
    input_eval(directory, csv_file) # check input
    individual_samples, ancestor_clone_pairs, anc_mult_clone_groups, set_of_pools, acp_trios, pipeline = read_csv_file(csv_file)

    # empty sets or lists initialized above evaluate to False; if any have members created by previous function, "if" statement below is True and dependent function is executed
    pileups_not_found = check_pileups(individual_samples)

    if not pileups_not_found: # pileup files need to be created
        print "All necessary pileup files were found. Run script 2 to continue analysis."
        sys.exit(0)
    else: fastq_single, fastq_pairs = find_fastq_files(pileups_not_found)
   
    if fastq_single:
        pileup_status_single = fastq_to_pileup(fastq_single)

    if fastq_pairs: 
        pileup_status_pairs = fastq_to_pileup(fastq_pairs) # make pileups

    if pileup_status_single or pileup_status_pairs != 0: # error while making pileup files; status is 0 if function not performed
        print "Error encountered while executing pileup function. Check sbatch output files."
        sys.exit(1)
    else:
        print "Run script 2 to generate varscan files when all pileups have been made."
        print ("To receive an email when pileups have been completed, confirm that your email address is in line 8 of a copy of script email_notification.sh. Then type the following into the command line, substituting in the sbatch job ID and full path to script:") 
        print "sbatch --dependency=afterok:<jobid> </your/path/to/email_notification.sh>"
        print "To check job status, type squeue -u <username> -j<jobid> into the command line."
        sys.exit(0)


if __name__ == "__main__": 
    main()

        
        
        
        
        
