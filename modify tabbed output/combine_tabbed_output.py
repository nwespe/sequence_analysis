import fnmatch, os, sys, string

''' input is from "run_combine_output.py - receive a list of tabbed_output files and output file name
'''
## modified 3/15/16 for comparing sawc files of various clones w/ same pop
## modified 3/18/16 for cmc file input

def main(arg0, arg1, arg2):  # mod - added arg0
    #separate gene from mutation files in list
    pop = arg0  # mod
    gene_files = {}
    mutation_files = {}
    file_dict = arg1
    samples = sorted(file_dict.iterkeys())
    print pop, samples
    for sample, files in file_dict.iteritems():
        for item in files:
            if fnmatch.fnmatch(item, '*gene.txt'):
                gene_files[sample] = item
            elif fnmatch.fnmatch(item, '*mutation.txt'):
                mutation_files[sample] = item

    mutation_dict = {}
    gene_dict = {}
    file_num = 0
    for sample in samples:
        file_num += 1
        input_file = open(gene_files[sample],'r')
        for line in input_file:
            entries = line.split('\t')
            keyname = tuple((entries[2], entries[3])) # chr_num, position; first row is keyed by 'chr_num', 'position'
            gene_key = tuple((entries[0], entries[2], entries[3])) # gene_name, chr_num, position; first row is keyed by 'gene_name', 'chr_num', 'position'
            if keyname in mutation_dict:
                # continue
                mutation_dict[keyname][3][sample] = entries[8]  # mod - now add sample name and fraction for output
            else:
                # mutation_dict[keyname] = entries[5:10]  # get list of entries -> 5 items
                mutation_dict[keyname] = entries[5:8]  # mod from [5:10] to [5:8] (now 3 items)
                fraction_dict = {sample: entries[8]}  # mod - now add sample name and fraction for output
                mutation_dict[keyname].append(fraction_dict)  # mod
            if gene_key in gene_dict:
                old_mut_type = gene_dict[gene_key]
                new_mut_type = entries[10].rstrip('\n')
                if new_mut_type == old_mut_type: continue
                elif new_mut_type == 'nonsegregating': continue
                elif old_mut_type == 'nonsegregating':
                    gene_dict[gene_key] = new_mut_type
            else:
                gene_dict[gene_key] = entries[10].rstrip('\n')

        #print mutation_dict['4','1163347']

    first_row = ('chr_num', 'position')
    column_list = []
    for sample in samples:  # in mod, samples is list of clones - only info from mutation file is percent in pop
        column_list.append(sample)
        read_file = open(mutation_files[sample], 'r')
        for line in read_file:
            entries = line.split('\t')
            keyname = tuple((entries[0], entries[1]))
            if keyname == first_row:
                # mutation_dict[first_row].append(sample)
                continue
            if keyname in mutation_dict:
            #    mutation_dict[keyname].append((sample, entries[3])) # add tuple of (pool name, seg percent) to list
                if pop in mutation_dict[keyname][3]:  # mod
                    continue  # mod
                else:   # mod
                    mutation_dict[keyname][3][pop] = entries[3]  # mod - add pop, seg percent to fraction_dict

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
    if os.path.isfile(output_name): os.remove(output_name)
    outputfile = open(output_name,'w')

    # outputfile.write('\t'.join(first_row) + '\t' + '\t'.join(mutation_dict[first_row]))
    outputfile.write('\t'.join(first_row) +'\t'+ '\t'.join(mutation_dict[first_row][0:3]))  # mod
    column_list.sort()
    column_list.append(pop)  # mod
    outputfile.write('\t'+ '\t'.join(column_list))  # mod
    i = 0
    while i < 4: # max_genes:
        i += 1
        outputfile.write('\tgene_name_' + str(i) + '\tmutation_type_' + str(i))
    outputfile.write('\n')
    mutation_dict.pop(first_row)

    # use column_list to match seg_percent with correct sample -> tuples start w entry 5 in mutation_dict list
    for keyname, entry in mutation_dict.iteritems():
        # outputfile.write('\t'.join(keyname) + '\t' + '\t'.join(entry[0:5]))
        outputfile.write('\t'.join(keyname) + '\t' + '\t'.join(entry[0:3]))  # mod
        # percent_dict = {} #initialize / reblank dictionary
        # for sample, percent in entry[5:]:
        #    percent_dict[sample] = percent # fill dictionary with sample percents for current gene entry
        for column in column_list:
            # if column in percent_dict: outputfile.write('\t' + percent_dict[column])
            if column in entry[3]:
                outputfile.write('\t' + entry[3][column])  # mod
            else:
                outputfile.write('\t')
        if keyname in position_dict:
            for (gene_name, mutation_type) in position_dict[keyname]:
                outputfile.write('\t' + gene_name + '\t' + mutation_type)
        outputfile.write('\n')

if __name__=='__main__':
    sys.exit(main(arg1, arg2))