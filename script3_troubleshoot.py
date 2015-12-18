#!/usr/bin/env python
# -*- coding: utf-8 -*-
# script3_troubleshoot.py

__author__ = 'nwespe'

#WTF

import csv
import re
import fnmatch
import os

directory = sys.argv[1]
csv_file = sys.argv[2]
csv_open = open(csv_file, 'rU')
csv_reader = csv.reader(csv_open)
row1 = next(csv_reader)

pool_regex = re.compile('(Pool).*')
pools = ([p.group(0) for cell in row1 for p in [pool_regex.search(cell)] if p])
num_pools = len(pools)
print num_pools

snp_files = []
for r, d, f in os.walk(directory):
    for item in f:
        if fnmatch.fnmatch(item, '*.snp'):
            snp_files.append(os.path.abspath(os.path.join(r, item)))

for pool in pools:
    item_pattern = '*'+pool+'*'
    matched_snp = fnmatch.filter(snp_files, item_pattern)
    print matched_snp
    remove_snp = fnmatch.filter(matched_snp, '*_vs_*') # need to remove anc_vs_pool files if they exist
    print remove_snp
    snp_file = set(matched_snp) - set(remove_snp)
    print snp_file
    # matched_indel = fnmatch.filter(indel_files, item_pattern)
    # print matched_indel
    # remove_indel = fnmatch.filter(matched_indel, '*_vs_*')
    # print remove_indel # check issue
    # indel_file = set(matched_indel) - set(remove_indel)
    # print indel_file