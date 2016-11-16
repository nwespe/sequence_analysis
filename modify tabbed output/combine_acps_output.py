import os
import sys
import fnmatch

__author__ = 'nmcollin'


class Analysis:
    """ ancestor-clone-pool-condition mutant analysis
    """

    def __init__(self, name, gene_file, mutation_file):
        self.name = name
        self.gene_file = gene_file
        self.mutation_file = mutation_file
        self.ancestor, self.clone, self.pool = self.name.split('_', 2)
        self.pool_strain, self.condition = self.pool.split('_', 1)

    def fill_mutation_dict(self, mutation_dict):
        for line in open(self.mutation_file, 'r'):
            entries = line.split('\t')
            if entries[0] == 'chr_num':
                continue  # first line
            chr_num = int(entries[0])
            position = int(entries[1])
            snp_indel = entries[2]
            seg_percent = entries[3]
            mutation_dict[(chr_num, position)] = Mutation(chr_num, position, snp_indel, self, seg_percent)
        return mutation_dict

    def parse_mutation_file(self, mutation_dict):
        for line in open(self.mutation_file, 'r'):
            entries = line.split('\t')
            if entries[0] == 'chr_num':
                continue  # first line
            chr_num = int(entries[0])
            position = int(entries[1])
            snp_indel = entries[2]
            seg_percent = int(entries[3])
            if (chr_num, position) in mutation_dict:
                mutation = mutation_dict[(chr_num, position)]  # call Mutation object
                mutation.seg_dict[self.condition] = seg_percent
            else:
                mutation_dict[(chr_num, position)] = Mutation(chr_num, position, snp_indel, self, seg_percent)
        return mutation_dict

    def parse_gene_file(self, mutation_dict, global_gene_dict):
        for line in open(self.gene_file, 'r'):
            entries = line.split('\t')
            if entries[0] == 'gene_name':
                continue  # first line
            chr_num = int(entries[2])
            position = int(entries[3])
            ref_base = entries[6]
            read_base = entries[7]
            clone_fraction = int(entries[8])
            mutation_type = entries[10].rstrip('\n')
            gene_name = entries[0]
            mutation = mutation_dict[(chr_num, position)]
            mutation.add_info_from_gene_file(ref_base, read_base, gene_name, mutation_type, clone_fraction)
            if gene_name in global_gene_dict:
                gene = global_gene_dict[gene_name]  # call Gene object
                gene.add_mutation(mutation)  # either create new entry in gene_mutation_dict or add Mutation to list
            else:
                global_gene_dict[gene_name] = Gene(gene_name, chr_num, mutation)
                # initialize Gene object with current Mutation
        return mutation_dict, global_gene_dict


class Comparison:
    """ collection of analyses using same ancestor/clone/pool trio
    """

    def __init__(self, a, c, p, condition, analysis):
        self.ancestor = a
        self.clone = c
        self.pool = p
        self.conditions = [condition]
        self.analyses = [analysis]
        self.mutation_dict = {}

    def collect_mutation_info(self, global_gene_dict):
        self.mutation_dict = self.analyses[0].fill_mutation_dict(self.mutation_dict)
        # fill mutation_dict from one mutation_file
        for analysis in self.analyses[1:]:
            self.mutation_dict = analysis.parse_mutation_file(self.mutation_dict)
            # get mutations/seg info from each additional Analysis mutation_file
        for analysis in self.analyses:
            self.mutation_dict, global_gene_dict = analysis.parse_gene_file(self.mutation_dict, global_gene_dict)
        # get mutation info from one gene_file, also add to global_gene_dict
        return global_gene_dict

    def write_comp_output(self, directory):
        output_basename = str(self.ancestor) + '_' + str(self.clone) + '_' + str(self.pool) + '_comparisons.txt'
        output_directory = os.path.abspath(directory)
        output_filename = os.path.abspath(os.path.join(output_directory, output_basename))
        sample_list = [self.clone] + self.conditions
        with open(output_filename, 'w') as output:
            output.write('chr_num \t position \t' + '\t'.join(sample_list) + '\t snp_indel \t ref_base \t read \t '
                         'gene_1 \t mut_type_1 \t gene_2 \t mut_type_2 \t gene_3 \t mut_type_3 \n')
            for (chr_num, position), mutation in sorted(self.mutation_dict.iteritems(), key = lambda x: x[0]):
                output.write(str(chr_num) + '\t' + str(position) + '\t')
                for sample in sample_list:
                    output.write(str(mutation.seg_dict[sample]) + '\t')
                output.write(mutation.snp_indel + '\t' + mutation.ref + '\t' + mutation.read)
                for gene, mut_type in mutation.gene_dict.iteritems():
                    output.write('\t' + gene + '\t' + mut_type)
                output.write('\n')


class Mutation:
    """
    one instance of Mutation class created by the Analysis class for each line in mutation_file of an analysis
    """

    def __init__(self, chr_num, position, snp_indel, analysis, seg_percent):
        self.seg_dict = {}
        self.gene_dict = {}  # gene_dict, ref, read and clone_fraction will only exist if mutation hits gene feature
        self.ref = 'unknown'
        self.read = 'unknown'
        self.chr_num = chr_num
        self.position = position
        self.snp_indel = snp_indel
        self.mut_sum = self.snp_indel + ': ' + self.ref + ' to ' + self.read
        self.clone, self.pool = analysis.clone, analysis.pool_strain
        self.seg_dict[analysis.condition] = seg_percent
        self.seg_dict[self.clone] = 'unknown'

    def add_info_from_gene_file(self, ref, read, gene, mutation_type, clone_fraction):
        if self.ref == 'unknown':
            self.ref = ref
            self.read = read
            self.mut_sum = self.snp_indel + ': ' + self.ref + ' to ' + self.read
            self.gene_dict[gene] = mutation_type
            self.seg_dict[self.clone] = clone_fraction
        elif gene in self.gene_dict and self.gene_dict[gene] != 'nonsegregating':
            return
        else:
            self.gene_dict[gene] = mutation_type


class Gene:
    """ called by parse_gene_file to create and add to global_gene_dict for summary output
        issue: only want to add strain info for clones, not pools
        currently have one Mutation object per analysis, not per clone
    """

    def __init__(self, name, chr_num, mutation):
        self.name = name
        self.gene_mutation_dict = {}
        self.hits = 1
        self.num_mutations = 1
        self.max_mut_hits = 1
        self.chr_num = chr_num
        self.gene_mutation_dict[(mutation.position, mutation.ref, mutation.read)] = \
            [mutation.mut_sum, mutation.gene_dict[self.name], (mutation.clone, mutation.seg_dict[mutation.clone])]
        # add mut_sum, mutation type and first mutation - only add clone name and fraction value

    def add_mutation(self, mutation):
        key = (mutation.position, mutation.ref, mutation.read)
        if key in self.gene_mutation_dict:
            if (mutation.clone, mutation.seg_dict[mutation.clone]) not in self.gene_mutation_dict[key]:
                self.gene_mutation_dict[key].append((mutation.clone, mutation.seg_dict[mutation.clone]))
                self.hits += 1
            # check mutation type and update in gene_mutation_dict if not nonsegregating
            if mutation.gene_dict[self.name] != 'nonsegregating':
                self.gene_mutation_dict[key][1] = mutation.gene_dict[self.name]
        else:
            self.gene_mutation_dict[key] = [mutation.mut_sum, mutation.gene_dict[self.name],
                                            (mutation.clone, mutation.seg_dict[mutation.clone])]
            self.num_mutations += 1
            self.hits += 1

    def get_max_mut_hits(self):
        position_hits = []
        for key, value in self.gene_mutation_dict.iteritems():
            position_hits.append(len(value))
        self.max_mut_hits = max(position_hits)
        return self.max_mut_hits


def get_analysis_files(directory):  # create list of Analysis objects
    analyses = []
    for item in os.listdir(directory):
        currdir = os.path.join(directory, item)
        if os.path.isdir(currdir):
            name = item
            gene_file = False
            mutation_file = False
            for subitem in os.listdir(currdir):
                if fnmatch.fnmatch(subitem, '*gene.txt'):
                    gene_file = os.path.join(currdir, subitem)
                elif fnmatch.fnmatch(subitem, '*mutation.txt'):
                    mutation_file = os.path.join(currdir, subitem)
            if gene_file and mutation_file:
                analyses.append(Analysis(name, gene_file, mutation_file))

    return analyses  # each analysis folder is in list, is an instance of class Analysis and has attributes


def create_comparison_dict(analyses):  # create output file with data from all 3 analyses for a given pool strain
    comparison_dict = {}  # (ancestor, clone, pool_strain) = Comparison object
    for analysis in analyses:
        a = analysis.ancestor
        c = analysis.clone
        p = analysis.pool_strain
        if (a, c, p) in comparison_dict:  # if Comparison object already exists for acp trio
            group = comparison_dict[(a, c, p)]  # call Comparison object
            group.conditions.append(analysis.condition)  # add condition to Comparison.conditions
            group.analyses.append(analysis)  # add Analysis object to list of analyses in Comparison
        else:
            comparison_dict[(a, c, p)] = Comparison(a, c, p, analysis.condition, analysis)  # create Comparison object

    return comparison_dict


def create_output_files(comparison_dict, global_gene_dict, current_dir):
    # parse input files and write output file for each comparison
    for acp_trio, group in comparison_dict.iteritems():
        global_gene_dict = group.collect_mutation_info(global_gene_dict)
        group.write_comp_output(current_dir)

    # write output file for global gene dict

    with open('gene_summary.txt', 'w') as output:
        output.write('hits_in_gene \t num_samples_w_mutation \t gene \t chr_num \t position \t mutation \t mut_type')

        # need to add enough column headings for max # of samples with same mutation
        # (i.e. # of entries in gene_mutation_dict for a given (position, read)
        gene_mut_hit_values = []
        for name, gene in global_gene_dict.iteritems():
            gene_mut_hit_values.append(gene.get_max_mut_hits())
        max_samples = max(gene_mut_hit_values)
        i = 1
        while i <= max_samples:
            output.write('\t sample_' + str(i) + '\t fraction_' + str(i))
            i += 1
        output.write('\n')

        # sort by hits attribute of Gene object (value)
        for (name, gene) in sorted(global_gene_dict.iteritems(), key = lambda x: x[1].hits, reverse=True):
            for key, mut_sample_list in sorted(gene.gene_mutation_dict.iteritems(), key = lambda x: x[0][0]):
                output.write(
                    str(gene.hits) + '\t' + str(len(mut_sample_list)-2) + '\t' + name + '\t' + str(gene.chr_num) +
                    '\t' + str(key[0]) + '\t' + mut_sample_list[0] + '\t' + mut_sample_list[1])
                # key[0] = position, last two entries are mut_sum and mutation_type
                for (strain, fraction) in sorted(mut_sample_list[2:], key = lambda x: x[0]):
                    output.write('\t' + strain + '\t' + str(fraction))
                output.write('\n')


def main():
    global_gene_dict = {}
    directory = os.path.abspath(sys.argv[1])
    analyses = get_analysis_files(directory)
    comparison_dict = create_comparison_dict(analyses)
    create_output_files(comparison_dict, global_gene_dict, directory)


if __name__ == "__main__":
    main()
