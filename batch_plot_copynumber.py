#!/usr/bin/env python
# -*- coding: utf-8 -*-
# run_plot_copynumber.py

# This script will execute the plot_copynumber.py script on the Harvard Odyssey cluster for samples specified by an
# input csv file.

import os, sys, re, fnmatch, csv, glob, argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

script_directory = os.path.dirname(sys.argv[0])

parser = argparse.ArgumentParser(description="Specify a directory where .copynumber files are present and a csv file "
                                             "describing samples to be analyzed.")
parser.add_argument('-i', '--input_directory', required=True,
                    help='Full path of directory containing .copynumber files. Files can be in subdirectories of this '
                         'directory.')
parser.add_argument('-f', '--csv_file', required=True, help='Full path to csv file containing names of samples.')
parser.add_argument('-o', '--output_directory', required=False, default='same',
                    help='Full path of directory for output files.')
parser.add_argument('-w', '--window_size', required=False, default=1000, help='Window size for smoothing of coverage.')
args = parser.parse_args()

input_dir = args.input_directory
output_dir = args.output_directory
if output_dir == 'same':
    output_dir = input_dir
csv_file = args.csv_file
window_size = args.window_size


def input_eval():

    """ Determine if inputs are a directory and a readable csv file
    """
    try:
        os.path.isdir(input_dir)
    except:
        print input_dir + " is not a directory"
        sys.exit(1)

    print "Now working on files in " + input_dir


def read_csv_file():

    """ Read csv file of sample names
    """

    csv_open = open(csv_file, 'rU')
    csv_reader = csv.reader(csv_open)
    row1 = next(csv_reader)

    clone_regex = re.compile('(Clone).*')
    clones = ([c.group(0) for cell in row1 for c in [clone_regex.search(cell)] if c])
    num_clones = len(clones)

    ancestor_index = row1.index('Ancestor')
    clone_indices = []
    for c in clones:
        clone_indices.append(row1.index(c))

    # make sets to pass to other functions
    individual_samples = set()  # set of all samples
    ancestor_clone_pairs = set()  # pairs to find copynumber files for
    for row in csv_reader:
        ancestor = row[ancestor_index]
        individual_samples.add(ancestor)
        list_of_clones = []
        for c in clone_indices:
            try:
                clone = row[c]
                if clone:
                    individual_samples.add(clone)
                    ancestor_clone_pairs.add((ancestor, clone))
                    list_of_clones.append(clone)
            except IndexError:  # no more clones in row
                break

    print "Received directions to analyze " + str(len(individual_samples)) + \
          " unique samples and " + str(len(ancestor_clone_pairs)) + " anc-clone pairs. "

    return individual_samples, ancestor_clone_pairs


def find_cn_files(pairs):

    """ Find copynumber files for ancestor_clone_pairs
    """
    # find all copynumber files in directory tree:
    cn_files = []
    for r, d, f in os.walk(input_dir):
        for item in f:
            if fnmatch.fnmatch(item, '*.copynumber'):
                cn_files.append(os.path.abspath(os.path.join(r, item)))

    # make dictionary of (key=anc_vs_clone, value=(ancestor, clone))
    # anc_vs_clone name is used to search copynumber files
    anc_clone_dict = {}
    for pair in pairs:  # ancestor_clone_pairs is a set of tuples
        anc_clone_name = str(pair[0] + '_vs_' + pair[1])  # name of analysis
        anc_clone_dict[anc_clone_name] = pair

    # create dictionary of (key = anc_vs_clone, value=(anc_vs_clone.copynumber))
    cn_dict = {}
    analyzed_pairs = set()
    missing_pairs = set()
    for key, pair in anc_clone_dict.iteritems():
        key_pattern = '*' + key + '*'
        matched_cn = fnmatch.filter(cn_files, key_pattern)
        if not matched_cn:
            missing_pairs.add(key)  # file not found
        else:
            analyzed_pairs.add(pair)
            cn_dict[key] = matched_cn[0]  # filename sent to plot_copynumber.py

    return missing_pairs, analyzed_pairs, cn_dict


def plot_copynumber(input_file, output_dir):

    filename = os.path.basename(input_file)
    print 'Now plotting ' + str(filename)
    name = os.path.splitext(filename)[0]  # get the names of the two genomes compared
    anc_text, sep, evo_text = name.partition("_vs_")

    #look for anc and evo individual plot files in output directory
    skip_anc, skip_evo = False, False
    if glob.glob(os.path.join(output_dir, anc_text + '_copynumber.pdf')): skip_anc = True
    if glob.glob(os.path.join(output_dir, evo_text + '_copynumber.pdf')): skip_evo = True

    chrom_dict = {}
    with open(input_file) as datafile:
        # The first line of the file is the header line. Store in variable `header,' which we will ignore.
        header = datafile.readline()
        # Read the file one line at a time
        for line in datafile:
            chrom = line.split()[0]  # The first position in the line gives the chromosome
            # If we don't already have an entry for the chromosome in our dictionary, create an entry.
            # The entry consists of four empty lists
            if chrom not in chrom_dict:
                chrom_dict[chrom] = ([], [], [], [], [])
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

    # calculate the median of depth across the whole genome
    anc_genome_median = np.median(np.array(anc_genome_list))
    evo_genome_median = np.median(np.array(evo_genome_list))

    chrom_list = sorted(chrom_dict.keys())
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

        W = window_size  # Calculate smoothed coverage over given window size (default 1000)
        total_bin_size = np.zeros((start[-1] / W) + 1)  # initializes array for smoothed coverage values
        anc_smoothed = np.zeros_like(total_bin_size)  # initializes another array with same shape
        evo_smoothed = np.zeros_like(total_bin_size)
        diff_smoothed = np.zeros_like(total_bin_size)

        for i in range(len(normalized_anc_depth)):  # for each entry in norm_anc_depth array
            total_bin_size[start[i] / W] += bin_size[i]
            anc_smoothed[start[i] / W] += normalized_anc_depth[i] * bin_size[i]
            evo_smoothed[start[i] / W] += normalized_evo_depth[i] * bin_size[i]
            diff_smoothed[start[i] / W] += depth_diff[i] * bin_size[i]
        anc_smoothed /= total_bin_size  # this method assumes all bins are equal size
        evo_smoothed /= total_bin_size
        diff_smoothed /= total_bin_size

        x_vals.append(W * np.arange(len(anc_smoothed)))
        y_vals_anc.append(anc_smoothed)
        y_vals_evo.append(evo_smoothed)
        y_vals_diff.append(diff_smoothed)
        all_y_vals = {anc_text: y_vals_anc, evo_text: y_vals_evo, name: y_vals_diff}
        if skip_anc: all_y_vals.pop(anc_text)
        if skip_evo: all_y_vals.pop(evo_text)

    # centromere locations in list form:
    centromeres = [(151465 - W, 151582 + W), (238207 - W, 238323 + W), (114385 - W, 114501 + W),
                   (449711 - W, 449821 + W), (151987 - W, 152104 + W), (148510 - W, 148627 + W),
                   (496920 - W, 497038 + W), (105586 - W, 105703 + W), (355629 - W, 355745 + W),
                   (436307 - W, 436425 + W), (440129 - W, 440246 + W), (150828 - W, 150947 + W),
                   (268031 - W, 268149 + W), (628758 - W, 628875 + W), (326584 - W, 326702 + W),
                   (555957 - W, 556073 + W)]

    for key, y_vals in all_y_vals.iteritems():  # create pdf files of chromosome plots
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

        fig, ax = plt.subplots(4, 4, sharey='row')
        fig.text(0.5, 0.01, 'chromosome coordinate', ha='center', va='center', fontsize=8)
        fig.text(0.01, 0.5, 'normalized read depth', ha='center', va='center', rotation='vertical', fontsize=8)

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
        fig.tight_layout()
        pdf_pages.savefig(fig)
        pdf_pages.close()
        plt.close(fig)

def main(): # run analyses of copynumber files using functions defined above

    """ check input,
        extract data from csv file
        find copynumber files created by varscan
        run plot_copynumber.py script ALT: add plot_copynumber code to this script
    """
    # initialize lists
    individual_samples = set()
    ancestor_clone_pairs = set()
    missing_pairs = set()
    analyzed_pairs = set()

    input_eval()  # check input
    individual_samples, ancestor_clone_pairs = read_csv_file()
    missing_pairs, analyzed_pairs, cn_dict = find_cn_files(ancestor_clone_pairs)

    if missing_pairs:
        print 'Varscan copynumber analysis files missing for pairs '+', '.join(str(f) for f in missing_pairs)
        print "Run copynumber.py script to generate varscan files or script 1 to generate pileup files."
    else: print 'All varscan copynumber analysis files have been found.'

    for input_file in cn_dict.itervalues():  # cn_dict = {'Anc_vs_Evo': /Anc_vs_Evo.copynumber}
        plot_copynumber(input_file, output_dir)

    print 'Plots are complete.'
    sys.exit(0)

if __name__ == "__main__":
    main()
