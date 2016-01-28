#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sequence_analysis.py

""" Description: Create snp and indel files using VarScan.
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

parser = argparse.ArgumentParser(description = "Tell me where to find pileup files, where to put the varscan files, and what samples to compare.")
parser.add_argument('-i', '--pileup', required=True, help='Full path of directory containing pileup files. Files can be in subdirectories of this directory.')
parser.add_argument('-d', '--directory', default='same', help='Full path of directory for file output. If omitted, output directory is same as pileup directory.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing names of ancestors, clones, and/or pools. See example documents for formatting.')
args = parser.parse_args()

pileup_directory = args.pileup
directory = args.directory
if directory == 'same':
    directory = pileup_directory
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


def check_snp_indel(pairs, pools):
    
    """ Check for snp and indel files for ancestor_clone_pairs, set_of_pools
    """
    # find all snp and indel files in directory tree:
    indel_files = []
    snp_files = []
    for r, d, f in os.walk(directory):
        for item in f:
            if fnmatch.fnmatch(item, '*.indel'):
                indel_files.append(os.path.abspath(os.path.join(r, item)))
            if fnmatch.fnmatch(item, '*.snp'):
                snp_files.append(os.path.abspath(os.path.join(r, item)))
                
    if not indel_files and not snp_files:  # no files yet
        samples = set() # initialize set of individual samples to find pileups for
        for pair in pairs: 
            for sample in pair: samples.add(sample)
        samples.update(pools) # combine all individual samples into one set
        print 'Now finding pileup files for ' + str(len(samples)) + ' samples'
        return samples, pairs, pools

    # below script should only run if any snp or indel files were found
    
    # make dictionary of (key=anc_vs_clone, value=(ancestor, clone))
    # anc_vs_clone name is used to search snp and indel files
    anc_clone_dict = {}
    for pair in pairs: # ancestor_clone_pairs is a set of tuples
        anc_clone_name = str(pair[0] + '_vs_' + pair[1]) # name of analysis
        anc_clone_dict[anc_clone_name] = pair

    # create sets of pairs, pools not already analyzed
    filtered_pairs = set()
    analyzed_pairs = set()
    for key, pair in anc_clone_dict.iteritems():
        key_pattern = '*'+key+'*'
        matched_snp = fnmatch.filter(snp_files, key_pattern)
        matched_indel = fnmatch.filter(indel_files, key_pattern)
        if not matched_snp or not matched_indel: filtered_pairs.add(pair) # missing file(s)
        else: analyzed_pairs.add(pair) # both files exist

    filtered_pools = set()
    analyzed_pools = set()
    for pool in pools:
        item_pattern = '*'+pool+'*'
        matched_snp = fnmatch.filter(snp_files, item_pattern)
        remove_snp = fnmatch.filter(matched_snp, '_vs_')
        snp_file = set(matched_snp) - set(remove_snp)
        matched_indel = fnmatch.filter(indel_files, item_pattern)
        remove_indel = fnmatch.filter(matched_indel, '_vs_')
        indel_file = set(matched_indel) - set(remove_indel)
        if not snp_file or not indel_file: filtered_pools.add(pool) # missing file(s)
        else: analyzed_pools.add(pool) # both files exist

    filtered_samples = set() # initialize set of individual samples to find pileups for
    for pair in filtered_pairs: 
        for sample in pair: filtered_samples.add(sample)
    filtered_samples.update(filtered_pools) # combine all individual samples into one set

    return filtered_samples, filtered_pairs, filtered_pools


def find_pileups(samples, pairs, pools):

    """ Get sets of ancestor-clone pairs and pools and find pileup files in filepath.
    """
    
    pileup_files = []
    for r, d, f in os.walk(pileup_directory):
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
# output directory not being used by varscan function?
 
    pool_pileups = []
    for pool in pools:
        if pool in pileups_dict:
            pool_pileup = pileups_dict[pool]
            output_dir = os.path.abspath(os.path.join(directory, 'varscan_files', pool))
            if not os.path.exists(output_dir): os.makedirs(output_dir)
            pool_pileups.append((pool_pileup, output_dir))

    if pileups_not_found: print 'Pileups not found for: ' + ', '.join(str(f) for f in pileups_not_found)
    if pileups_not_unique: print 'ID names not unique for: ' + ', '.join(str(f) for f in pileups_not_unique)

    return paired_pileups, pool_pileups, pileups_not_found #, pileups_not_unique


def batch_varscan_clones(paired_pileups): # paired_pileups is a list of tuples (anc.pileup, clone.pileup, anc_clone_name, output_dir)
    
    """ Run varscan somatic routine for ancestor-clone pairs from find_pileups
        function.
    """
    array_num = len(paired_pileups)-1
    bash_file = os.path.join(script_directory, 'array_scripts', 'varscan_clones.sh')
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    output = os.path.join(diag_dir, 'varscan_clones_%A-%a.out')
    error = os.path.join(diag_dir, 'varscan_clones_%A-%a.err')
    print "Now creating varscan files for " + str(len(paired_pileups)) + " ancestor-clone pairs"
    arg_list = []
    for pair in paired_pileups:
        arg_list.append(' '.join(pair))
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file, script_directory]  
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch varscan_clones.sh executed, return code is " + str(returncode)
    return returncode


def batch_varscan_pools(pool_pileups): # pool_pileups is a list of tuples (pool_pileup, output_dir)

    """ Run varscan pileuptosnp and pileuptoindel routines for pools from
        find_pileups function.
    """
    array_num = len(pool_pileups)-1
    bash_file = os.path.join(script_directory, 'array_scripts', 'varscan_pools.sh')
    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    output = os.path.join(diag_dir, 'varscan_pools_%A-%a.out')
    error = os.path.join(diag_dir, 'varscan_pools_%A-%a.err')
    print "Now creating varscan files for " + str(len(pool_pileups)) + " pools"
    arg_list = []
    for pool in pool_pileups:
        arg_list.append(' '.join(pool))
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file, script_directory]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch varscan_pools.sh executed, return code is " + str(returncode)
    return returncode


def main(): # run analyses of pileup files using functions defined above

    """ check input,
        extract data from csv file
        check if snp and indel files exist
        run varscan function if snp, indel files don't exist - in array format
    """
   
    # initialize lists
    samples = set(); pairs = set(); pools = set()
    fastq_pairs = []
    pileups_not_found = set()

    # set default exit statuses
    var_clone_status = 0; var_pool_status = 0; new_snp_indel = 0
    
    input_eval(directory, csv_file) # check input
    individual_samples, ancestor_clone_pairs, anc_mult_clone_groups, set_of_pools, acp_trios, pipeline = read_csv_file(csv_file)

    # empty sets or lists initialized above evaluate to False; if any have members created by previous function, "if" statement below is True and dependent function is executed
    samples, pairs, pools = check_snp_indel(ancestor_clone_pairs, set_of_pools)

    if not samples:
        print "All varscan analyses have been performed. Run script 3 to continue analysis."
        sys.exit(0)
    
    paired_pileups, pool_pileups, pileups_not_found = find_pileups(samples, pairs, pools)

    if pileups_not_found:
        print "Not all pileup files were found. Run script 1 to generate pileups."
        sys.exit(1) # pileup files need to be created
        
    # now all required pileup files should exist - create snp, indel files if needed
    if paired_pileups:
        new_snp_indel = 1
        var_clone_status = batch_varscan_clones(paired_pileups)
    if pool_pileups:
        new_snp_indel = 1
        var_pool_status = batch_varscan_pools(pool_pileups)
    if var_clone_status != 0 or var_pool_status !=0: # error while making varscan files; status is 0 if function not performed
        print "Error encountered while making snp, indel files. Check error log for varscan jobs."
        sys.exit(1)

    if new_snp_indel == 1: # if either varscan function was run
        print "Varscan analyses being run. Run script 3 to generate mutant comparison output when all files have been made. Use "+pipeline+" analysis pipeline for this csv file."
        print "To receive an email when varscan files have been completed, confirm that your email address is in line 8 of a copy of script email_notification.sh. Then type the following into the command line, substituting in the sbatch job ID(s) and full path to script:"
        print "sbatch --dependency=afterok:<jobid1>:<jobid2> </your/path/to/email_notification.sh>"
        print "To check job status, type squeue -u <username> -j<jobid> into the command line."
        sys.exit(0)
    else: 
        print "Error encountered while executing varscan function."
        sys.exit(1)


if __name__ == "__main__": 
    main()

        
