import os, sys

''' input is from "run_combine_cmc.py - receive a list of tabbed_output files and output file name
'''
## modified 3/15/16 for comparing sawc files of various clones w/ same pop## modified 3/18/16 for cmc file input

def main(arg0, arg1, arg2):
    name = arg0  # mod
    gene_file = arg1
    print name

    mutation_dict = {}
    gene_dict = {}
    sample_set = set()

    input_file = open(gene_file,'rU')
    for line in input_file:
        entries = line.split('\t')
        sample = entries[4]
        sample_set.add(sample)
        keyname = tuple((entries[2], entries[3])) # chr_num, position; first row is keyed by 'chr_num', 'position'
        gene_key = tuple((entries[0], entries[2], entries[3])) # gene_name, chr_num, position; first row is keyed by 'gene_name', 'chr_num', 'position'
        if keyname in mutation_dict:
            mutation_dict[keyname][3][sample] = (entries[8], entries[9])  # add sample name and fraction, tot_reads for output
        else:
            mutation_dict[keyname] = entries[5:8]
            fraction_dict = {sample: (entries[8], entries[9])}  # create fraction dict with entry for each sample
            mutation_dict[keyname].append(fraction_dict)
        if gene_key not in gene_dict:
            gene_dict[gene_key] = entries[10].rstrip('\n')
    print sample_set

    # compile position_dict from gene_dict
    position_dict = {}
    for (gene_name, chr_num, position) in gene_dict.iterkeys():
        if (chr_num, position) in position_dict:
            position_dict[(chr_num, position)].append((gene_name, gene_dict[(gene_name, chr_num, position)]))
        else:
            position_dict[(chr_num, position)] = [(gene_name, gene_dict[(gene_name, chr_num, position)])]
    # get length of longest element
    #max_genes = max([len(value) for key, value in position_dict.iteritems()])

    output_name = arg2
    if os.path.isfile(output_name):
        os.remove(output_name)
        print "deleted previous cmc_combined output file"
    outputfile = open(output_name,'w')

    first_row = ('chr_num', 'position')
    outputfile.write('\t'.join(first_row) +'\t'+ '\t'.join(mutation_dict[first_row][0:3]))
    sample_set.remove('sample_name')
    column_list = sorted(list(sample_set))
    print column_list
    for sample in column_list:
        outputfile.write('\t'+sample+'_fraction'+'\t'+sample+'_reads')
    i = 0
    while i < 4:  # max_genes:
        i += 1
        outputfile.write('\tgene_name_' + str(i) + '\tmutation_type_' + str(i))
    outputfile.write('\n')
    mutation_dict.pop(first_row)

    # use column_list to match fraction with correct sample column
    for keyname, entry in mutation_dict.iteritems():
        outputfile.write('\t'.join(keyname) + '\t' + '\t'.join(entry[0:3]))
        for column in column_list:
            if column in entry[3]:
                outputfile.write('\t' + '\t'.join(entry[3][column]))  # each entry is tuple of fraction, tot_reads
            else:
                outputfile.write('\t\t')
        if keyname in position_dict:
            for (gene_name, mutation_type) in position_dict[keyname]:
                outputfile.write('\t' + gene_name + '\t' + mutation_type)
        outputfile.write('\n')

if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))