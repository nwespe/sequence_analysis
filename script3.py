#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sequence_analysis.py

""" Description
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

parser = argparse.ArgumentParser(description = "Tell me where to find bam and varscan files, where to put the analysis files, and what samples to compare.")
parser.add_argument('-i', '--input_files', required=True, help='Full path of directory containing bam and varscan files. Files can be in subdirectories of this directory.')
parser.add_argument('-d', '--directory', required=True, help='Full path of directory for analysis output.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing desired comparisons and names of ancestors, clones, pools. See sample documents for formatting.')
parser.add_argument('-m', '--mutation_percent', type=int, default=90, help='Percent of reads required for mutation in clone to be called true. Default is 90; recommend 30 for diploids.')
parser.add_argument('-s', '--segregation_percent', type=int, default=70, help='Percent of reads required for mutation in pool to be called segregating. Default is 70.')
args = parser.parse_args()

input_dir = args.input_files
directory = args.directory
csv_file = args.csv_file
mut_percent = args.mutation_percent
seg_percent = args.segregation_percent

def input_eval():
    
    """ Determine if inputs are a directory and a readable csv file
    """
    try:
        os.path.isdir(directory)
    except:
        print directory + " is not a directory"
        sys.exit(1)
##TODO    try:
##        csv_open = open(csv_file, 'rU')
##        csv_reader = csv.reader(csv_open)
##        #some other function here to check file?
##    except:
##        print csv_file + " is not readable csv file."
##        return 1
        
    print "Now working on files in " + directory


def read_csv_file():
    
    """ Read csv file of analyses to be done and create input for find_pileups
        function. First row must be header row with "Ancestor, Clone, Pool_x"
        Returns four sets: individual samples, ancestor-clone pairs, pools, and trios.
    """
# problem: can't read csv file with uneven numbers of entries in rows
# hack: format all empty cells as text in excel before saving as csv
    
    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)
    row1 = next(csv_reader)

    # get number of clones and pools
    clone_regex = re.compile('(Clone).*')
    clones = ([c.group(0) for cell in row1 for c in [clone_regex.search(cell)] if c])
    num_clones = len(clones)
    pool_regex = re.compile('(Pool).*')
    pools = ([p.group(0) for cell in row1 for p in [pool_regex.search(cell)] if p])
    num_pools = len(pools)

    # get indices for ancestor, clone(s), pool(s)
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
    else: print "Please check csv file and format for either multi-clone comparisons or ancestor-clone-pool trio comparisons." 

    print "Received directions to analyze " + str(len(individual_samples)) + \
          " unique samples, " + str(len(ancestor_clone_pairs)) + " anc-clone pairs, " \
          + str(len(anc_mult_clone_groups)) + " multi-clone groups, " \
          + str(len(set_of_pools)) + " pools, and " + str(len(acp_trios)) + " trios."

    return individual_samples, ancestor_clone_pairs, anc_mult_clone_groups, set_of_pools, acp_trios, pipeline 


def find_snp_indel(pairs, pools):
    
    """ Search for snp and indel files for ancestor_clone_pairs, set_of_pools
    """
    # find all snp and indel files in directory tree:
    indel_files = []
    snp_files = []
    for r, d, f in os.walk(input_dir):
        for item in f:
            if fnmatch.fnmatch(item, '*.indel'):
                indel_files.append(os.path.abspath(os.path.join(r, item)))
            if fnmatch.fnmatch(item, '*.snp'):
                snp_files.append(os.path.abspath(os.path.join(r, item)))

    # make dictionary of (key=anc_vs_clone, value=(ancestor, clone))
    # anc_vs_clone name is used to search snp and indel files
    anc_clone_dict = {}
    for pair in pairs: # ancestor_clone_pairs is a set of tuples
        anc_clone_name = str(pair[0] + '_vs_' + pair[1]) # name of analysis
        anc_clone_dict[anc_clone_name] = pair
    #print 'anc_clone_dict = ', anc_clone_dict

    # create sets of pairs, pools already analyzed
    # and dictionaries of (key = anc_vs_clone or pool, value=(.snp, .indel))
    snp_indel_dict = {}
    analyzed_pairs = set()
    missing_pairs = set()
    for key, pair in anc_clone_dict.iteritems():
        key_pattern = '*'+key+'*'
        matched_snp = fnmatch.filter(snp_files, key_pattern)
        matched_indel = fnmatch.filter(indel_files, key_pattern)
        if not matched_snp or not matched_indel: missing_pairs.add(key) # one or both files not found
        else:
            analyzed_pairs.add(pair)
            snp_indel_dict[key] = (matched_snp[0], matched_indel[0]) # filenames sent to mutantanalysis

    analyzed_pools = set()
    missing_pools = set()
    for pool in pools:
        item_pattern = '*'+pool+'*'
        matched_snp = fnmatch.filter(snp_files, item_pattern)
        remove_snp = fnmatch.filter(matched_snp, '*_vs_*') # need to remove anc_vs_pool files if they exist
        snp_file = set(matched_snp) - set(remove_snp)
        matched_indel = fnmatch.filter(indel_files, item_pattern)
        remove_indel = fnmatch.filter(matched_indel, '*_vs_*')
        indel_file = set(matched_indel) - set(remove_indel)
        if not matched_snp or not matched_indel: missing_pools.add(pool) # one or both files not found
        else:
            analyzed_pools.add(pool)
            snp_indel_dict[pool] = (list(snp_file)[0], list(indel_file)[0]) # save filenames for analysis
    return missing_pairs, missing_pools, snp_indel_dict, analyzed_pairs, analyzed_pools


def mutant_analysis(pipeline, snp_indel_dict, pairs, mult_clone_groups, pools, trios):

    # dict: {(anc_vs_clone or pool=(.snp, .indel))}
    # pairs: set((anc, clone))
    # mult_clone_groups: set((anc, [list_of_clones])
    # pools: set(pool)
    # trios: set((anc, clone, pool))

    # check mult_clone_groups for those with all members in pairs sets, i.e. ready to analyze
    checked_mc_groups = []
    for anc, group in mult_clone_groups:
        all_good = False
        for clone in group:
            if (anc, clone) in pairs: all_good = True
        if all_good: checked_mc_groups.append((anc, group))

    # check trios for those with all three members in pairs&pools sets, i.e. ready to analyze
    checked_trios = []
    for trio in trios:
        if (trio[0], trio[1]) in pairs and trio[2] in pools: checked_trios.append(trio)

    # find bam files
    # TODO: add checks for sample bam files - make sure there is only 1 per sample
    bam_files = []
    for r, d, f in os.walk(input_dir):
        for item in f:
            if fnmatch.fnmatch(item, '*final.bam'):
                bam_files.append(os.path.abspath(os.path.join(r, item)))
    # print "Found " + str(len(bam_files)) + " bam files"

    # arguments for mutantanalysis.py
    ref_seq = os.path.join(script_directory, 'mutantanalysis', 'ref_seq', 's288c_sgd.fa')
    ref_feat = os.path.join(script_directory, 'mutantanalysis', 'ref_seq', 's288c_ref_annot.txt')
    ref_chrom = os.path.join(script_directory, 'mutantanalysis', 'ref_seq', 's288c_chr_names.txt')

    diag_dir = os.path.join(directory, 'diagnostic')
    if not os.path.exists(diag_dir): os.makedirs(diag_dir)
    output = os.path.join(diag_dir, 'mutantanalysis_%A-%a.out')
    error = os.path.join(diag_dir, 'mutantanalysis_%A-%a.err')

    analysis_dir = 'analysis_output-' + pipeline +'-m'+ str(mut_percent) +'-s'+ str(seg_percent)
    arg_list = []
    if pipeline == "sawc":
        array_num = len(checked_trios)-1
        for trio in checked_trios:
            anc_name = trio[0]
            clone_name = trio[1]
            pool_name = trio[2]
            key = anc_name+'_vs_'+clone_name 
            title = anc_name +'_'+ clone_name +'_'+ pool_name
            output_dir = os.path.abspath(os.path.join(directory, analysis_dir, title))
            if not os.path.exists(output_dir): os.makedirs(output_dir)
            anc_bam = fnmatch.filter(bam_files, '*'+anc_name+'*')[0]
            clone_bam = fnmatch.filter(bam_files, '*'+clone_name+'*')[0]
            clone_snp = snp_indel_dict[key][0]
            clone_indel = snp_indel_dict[key][1]
            pool_bam = fnmatch.filter(bam_files, '*'+pool_name+'*')[0]
            pool_snp = snp_indel_dict[pool_name][0]
            pool_indel = snp_indel_dict[pool_name][1]
            sample = (title, anc_name, anc_bam, output_dir, clone_name, clone_bam,
                      clone_snp, clone_indel, pool_name, pool_bam, pool_snp, pool_indel) 
                      # create tuple containing args for mutantanalysis       
            arg_list.append(' '.join(sample))  # convert tuple to space-delimited string and add to list
    elif pipeline == "cmc":    
        array_num = len(checked_mc_groups)-1
        for anc_name, mc_group in checked_mc_groups:
            anc_bam = fnmatch.filter(bam_files, '*'+anc_name+'*')[0]
            clone_list = []
            names = [anc_name]
            num_clones = len(mc_group)
            for clone_name in mc_group:
                clone_bam = fnmatch.filter(bam_files, '*'+clone_name+'*')[0]
                key = anc_name+'_vs_'+clone_name 
                clone_snp = snp_indel_dict[key][0]
                clone_indel = snp_indel_dict[key][1]
                clone_info = (clone_name, clone_bam, clone_snp, clone_indel)
                clone_list.append(' '.join(clone_info))
                names.append(clone_name)
            title = '_'.join(names)
            output_dir = os.path.abspath(os.path.join(directory, analysis_dir, title))
            if not os.path.exists(output_dir): os.makedirs(output_dir)
            sample = [title, anc_name, anc_bam, output_dir, str(num_clones)]
            sample.extend(clone_list) 
            # create list containing args for mutantanalysis
            arg_list.append(' '.join(sample))  # convert list to space-delimited string and add to list
    bash_file = os.path.join(script_directory, 'array_scripts', 'segregant_analysis.sh')
    subcall_list = ['sbatch', '--array=0-'+str(array_num), '--error='+str(error), '--output='+str(output), bash_file, '-d', script_directory, '-p', pipeline, '-r', ref_seq, '-c', ref_chrom, '-f', ref_feat, '-m', str(mut_percent), '-s', str(seg_percent)]
    subcall_list.extend(arg_list)
    returncode = subprocess.call(subcall_list)
    print "sbatch segregant_analysis.sh executed, return code is " + str(returncode)
    return returncode


def main(): # run analyses of pileup files using functions defined above

    """ check input,
        extract data from csv file
        find snp & indel files created by varscan
        run final mutantanalysis.py script
    """
    # initialize lists
    ancestor_clone_pairs = set(); anc_mult_clone_groups = []; 
    set_of_pools = set(); acp_trios = set()
    missing_pairs = set(); missing__pools = set()
    analyzed_pairs = set(); analyzed_pools = set()

    input_eval() # check input
    individual_samples, ancestor_clone_pairs, anc_mult_clone_groups, set_of_pools, acp_trios, pipeline = read_csv_file()
    missing_pairs, missing_pools, snp_indel_dict, analyzed_pairs, analyzed_pools = find_snp_indel(ancestor_clone_pairs, set_of_pools)

    if missing_pairs or missing_pools:
        print 'Varscan analysis files missing for pairs '+', '.join(str(f) for f in missing_pairs)+' and pools '+', '.join(str(f) for f in missing_pools)
        print "Run script 2 to generate varscan files or script 1 to generate pileup files."
        sys.exit(1)
    else: print 'All varscan analysis files have been found.'   # all required snp, indel files should exist - run mutantanalysis script
    mutant_analysis(pipeline, snp_indel_dict, analyzed_pairs, anc_mult_clone_groups, analyzed_pools, acp_trios) 
    print "Mutant comparison analysis is running. To check status, type squeue -u <username> -j<jobid> into the command line."
    print ("To receive an email when mutant comparison analyses have been completed, confirm that your email address is in line 8 of a copy of script email_notification.sh. Then type the following into the command line, substituting in the sbatch job ID and full path to script:")
    print "sbatch --dependency=afterok:<jobid1> </your/path/to/email_notification.sh>"
    sys.exit(0)


if __name__ == "__main__": 
    main()

