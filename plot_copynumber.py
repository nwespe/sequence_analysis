#!/usr/bin/env python
# -*- coding: utf-8 -*-
# plot_copynumber.py

"""
This program will ultimately create 3 pdf files: 2 of plots showing copy number normalized to the genome
(one for ancestor, one for evolved) and one of plots showing evo/anc copy number
Script originally written by Daniel Rice 2016, modified by Phoebe Hsieh, Marco Fumasoni, and Nichole Wespe 2016
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys, os
import seaborn as sns
import argparse

# Usage:
# python plot_copynumber.py -f DATAFILE.copynumber -d /output_directory/

parser = argparse.ArgumentParser(description='Specify a Varscan copynumber file to be plotted.')
parser.add_argument('-d', '--output_directory', required=False, help='Full path of output directory.')
parser.add_argument('-f', '--cn_file', required=True, help='Full path to copynumber file.')
args = parser.parse_args()

input = args.cn_file
output_dir = args.output_directory

# Read in the start position and log2_ratio for each window, saving each chromosome separately
# Initialize a dictionary to hold the data, arranged by chromosome
chrom_dict = {}
with open(input) as datafile:
    # The first line of the file is the header line. Store in variable `header,' which we will ignore.
    header = datafile.readline()
    # Read the file one line at a time
    for line in datafile:
        chrom = line.split()[0] # The first position in the line gives the chromosome
        # If we don't already have an entry for the chromosome in our dictionary, create an entry.
        # The entry consists of four empty lists
        if chrom not in chrom_dict:
            chrom_dict[chrom] = ([],[],[],[],[])
        # Store the start position, coverage_ratio, normal_depth, tumor_depth, bin_size in their respective lists
        chrom_dict[chrom][0].append(int(line.split()[1]))
        chrom_dict[chrom][1].append(float(line.split()[-2]))
        chrom_dict[chrom][2].append(float(line.split()[-4]))
        chrom_dict[chrom][3].append(float(line.split()[-3]))
        chrom_dict[chrom][4].append(int(line.split()[3]))

# convert chrom_dict list entries into arrays and create full genome list
anc_genome_list = []
evo_genome_list = []
for key, value in chrom_dict.iteritems():
    start = np.array(value[0])
    cov_ratio = np.array(value[1])
    normal_depth = np.array(value[2])
    tumor_depth = np.array(value[3])
    bin_size = np.array(value[4])
    chrom_dict[key] = (start, cov_ratio, normal_depth, tumor_depth, bin_size)
    anc_genome_list.extend(value[2])
    evo_genome_list.extend(value[3])

#calculate the median of depth across the whole genome
anc_genome_median = np.median(np.array(anc_genome_list))
evo_genome_median = np.median(np.array(evo_genome_list))

chrom_list = sorted(chrom_dict.keys())
ch = 1
x_vals = []
y_vals_anc = []
y_vals_evo = []
y_vals_diff = []

for chrom in chrom_list:  # calculate smoothed coverage for all 3 plots
    start = chrom_dict[chrom][0]
    cov_ratio = chrom_dict[chrom][1]
    anc_depth = chrom_dict[chrom][2]
    evo_depth = chrom_dict[chrom][3]
    bin_size = chrom_dict[chrom][4]
    # Normalize by median of genome depth to control for difference sequencing depths
    normalized_anc_depth = anc_depth / anc_genome_median
    normalized_evo_depth = evo_depth / evo_genome_median
    depth_diff = (normalized_evo_depth) - (normalized_anc_depth)

    W = 1000  # Calculate smoothed coverage over window size of 1000
    total_bin_size = np.zeros((start[-1] / W) + 1)  # initializes array for smoothed coverage values, e.g. chr1 = 231
    anc_smoothed = np.zeros_like(total_bin_size)  # initializes another array with same shape
    evo_smoothed = np.zeros_like(total_bin_size)
    diff_smoothed = np.zeros_like(total_bin_size)

    for i in range(len(normalized_anc_depth)):  # for each entry in norm_anc_depth array
        total_bin_size[start[i]/W] += bin_size[i]
        anc_smoothed[start[i]/W] += normalized_anc_depth[i]*bin_size[i]
        evo_smoothed[start[i]/W] += normalized_evo_depth[i]*bin_size[i]
        diff_smoothed[start[i]/W] += depth_diff[i]*bin_size[i]
    anc_smoothed /= total_bin_size  # this method assumes all bins are equal size
    evo_smoothed /= total_bin_size
    diff_smoothed /= total_bin_size

    x_vals.append(W*np.arange(len(anc_smoothed)))
    y_vals_anc.append(anc_smoothed)
    y_vals_evo.append(evo_smoothed)
    y_vals_diff.append(diff_smoothed)
    name = os.path.splitext(os.path.basename(input))[0]  # get the names of the two genomes compared
    anc_text, sep, evo_text = name.partition("_vs_")
    all_y_vals = {anc_text: y_vals_anc, evo_text: y_vals_evo, name: y_vals_diff}

# centromere locations in list form:
centromeres = [(151465 - W, 151582 + W), (238207 - W, 238323 + W), (114385 - W, 114501 + W), (449711 - W, 449821 + W),
               (151987 - W, 152104 + W), (148510 - W, 148627 + W), (496920 - W, 497038 + W), (105586 - W, 105703 + W),
               (355629 - W, 355745 + W), (436307 - W, 436425 + W), (440129 - W, 440246 + W), (150828 - W, 150947 + W),
               (268031 - W, 268149 + W), (628758 - W, 628875 + W), (326584 - W, 326702 + W), (555957 - W, 556073 + W)]


for key, y_vals in all_y_vals.iteritems(): # create pdf files of chromosome plots
    if key == name:
        ymin = -3
        ymax = 3
        color_min = -1
        color_max = 1
    else:
        ymin = -2
        ymax = 4
        color_min = 0
        color_max = 2
    fig_fn = os.path.join(output_dir, key + '_copynumber.pdf')
    pdf_pages = PdfPages(fig_fn)

    fig, ax = plt.subplots(4,4, sharey ='row')
    fig.text(0.5, 0.01, 'chromosome coordinate', ha='center', va='center',fontsize=8)
    fig.text(0.01, 0.5, 'normalized read depth', ha='center', va='center', rotation='vertical',fontsize=8)

    for i, a in enumerate(ax.flatten()):
        a.vlines(centromeres[i][0], ymin, ymax, alpha=0.1)
        a.vlines(centromeres[i][1], ymin, ymax, alpha=0.1)
        a.plot(x_vals[i], y_vals[i], alpha=0.3)
        cp = y_vals[i]
        a.scatter(x_vals[i], y_vals[i], marker='.', s=10, linewidth='0', c=cp, vmin=color_min, vmax=color_max,
                  cmap='rainbow', zorder=10)
        a.set_ylim([ymin, ymax])
        a.yaxis.set_ticks(np.arange(ymin, (ymax + 0.5), 0.5))
        xmax = max(x_vals[i])
        a.set_xlim(-15000, xmax + 15000)
        a.xaxis.set_tick_params(labelsize=3)
        a.yaxis.set_tick_params(labelsize=3)
        a.set_title('chromosome ' + str(i + 1), fontsize=8)
    fig.set_tight_layout(True)
    pdf_pages.savefig(fig)
    pdf_pages.close()
    plt.close(fig)
