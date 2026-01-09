#!/usr/bin/env python3
# todoc
"""
biom2biotaviz
-----------

.. module:: biom2biotaviz
  :synopsis: Convert biom file to BiotaViz-style txt file

Script generates a BiotaViz-style tab-delimited txt file from a biom file (v1).
Output is written to <infile>.txt and <infile>.biotaviz.txt

Typical run::

    biom2biotaviz.py -i some_biom_file.biom1 -o BiotaViz.txt
    
Changes:
28/06/2020: String formating update, added 'd' to label_replace dictionary

Author: Jos Boekhorst
"""
# Import required functions
import sys
import os
from argparse import ArgumentParser


def usage():
    sys.stderr.write(f"Use: {sys.argv[0]} <infile>\n")


def read_txt(filename):
    infile = open(filename, 'r')
    text = infile.read().rstrip('\n')
    infile.close()
    return text


def remove_empty_terminals(taxon):
    """remove terminal taxonomy bits that are empty"""
    new_taxon = []
    taxon_split = taxon.split('; ')
    for element in taxon_split:
        if element.split('__')[-1] != "":
            new_taxon.append(element)
        else:
            break
    return '; '.join(new_taxon)


def read_OTU_table(filename):
    # first line is "from biom" comment, second line is header
    taxonomy = {}
    counts = {}
    collapsed = {}
    lines = read_txt(filename)
    print(lines)
    lines = lines.split('\n')
    samples = lines[1].split('\t')[1:-1]
    for line in lines[2:]:
        lineg = line.split('\t')
        OTU = lineg[0]
        if lineg[-1] == 'Unassigned':  # happens in not-quite-filtered-enough qiime2 data, for example the test case
            lineg[-1] = 'd__Unassigned'
        taxon = 'r__Root; ' + lineg[-1]
        taxon = taxon.replace("NA;", "k__;")  # NG_Tax weirdness
        taxon = remove_empty_terminals(taxon)
        counts[OTU] = {}
        taxonomy[OTU] = taxon
        if taxon not in collapsed:
            collapsed[taxon] = {}
        for i, sample in enumerate(samples):
            count = float(lineg[i + 1])
            counts[OTU][sample] = count
            if sample not in collapsed[taxon]:
                collapsed[taxon][sample] = 0
            collapsed[taxon][sample] += count
    return collapsed, taxonomy, samples


def get_new_number(taxon_to_trace, trace):
    """get a new number for a taxon trace"""
    parent = trace[:-1]
    # get current highest child
    current_children = [0]
    for element in taxon_to_trace.keys():
        element_parent = element[:-1]
        if element_parent == parent:
            current_children.append(int(taxon_to_trace[element].rstrip('.').split('.')[-1]))
    return max(current_children) + 1


def traces_from_taxonomy(collapsed):
    taxon_to_trace = {tuple(['r__Root']): '1.'}
    for taxon in collapsed:
        taxon_gs = taxon.split('; ')
        full_name = []
        for taxon_g in taxon_gs:
            label = taxon_g
            full_name.append(label)

        # Now get a code. Parents need to be inferred.
        for i in range(1, len(full_name)):
            tmp = full_name[:i + 1]
            parent = full_name[:i]
            if not tuple(tmp) in taxon_to_trace:
                new_nr = get_new_number(taxon_to_trace, tuple(tmp))
                new_trace = f"{taxon_to_trace[tuple(parent)]}{str(new_nr)}."
                taxon_to_trace[tuple(tmp)] = new_trace
    return taxon_to_trace


def inverse_dict(input_dict):
    # generate new dictionary by swapping key and value
    new_dict = {}
    for element in input_dict.keys():
        element_value = input_dict[element]
        new_dict[element_value] = element
    return new_dict


def infer_internal_counts(collapsed):
    """data was collapsed per taxon, now infer parent counts"""
    new_counts = {}
    for taxon in collapsed.keys():
        taxon_g = taxon.split('; ')
        for i in range(len(taxon_g)):
            taxon_sub = tuple(taxon_g[:i + 1])
            if taxon_sub not in new_counts.keys():
                new_counts[taxon_sub] = {}
            for sample in collapsed[taxon].keys():
                if sample not in new_counts[taxon_sub].keys():
                    new_counts[taxon_sub][sample] = 0
                new_counts[taxon_sub][sample] += collapsed[taxon][sample]
    return new_counts


# settings
description = 'Converts biom file or biom-style OTU table to biotaviz file. Output is printed to std, redirect to file with "biom convert -i infile.biom  > BiotaViz.txt". You may also be interested in JOS_clean_biom_txt.py.'

label_replace = {'r': 'no',
                 'k': 'domain',
                 'd': 'domain',
                 'p': 'phylum',
                 'c': 'class',
                 'o': 'order',
                 'f': 'family',
                 'g': 'genus',
                 's': 'species',
                 'sh': 'specieshypothesis',
                 't': 'variant'}

if __name__ == "__main__":

    parser = ArgumentParser(description=description, add_help=True)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-o', dest='outfile', help='name of output file', default="stdout")
    parser.add_argument('-t', dest='isText', action='store_true', help="input is OTU table, not biom")
    options = vars(parser.parse_args())
    infile = options['infile']

    if options['isText'] is not True:
        outfile = infile + '.txt'
        if os.path.isfile(outfile):
            sys.stderr.write(
                "Outfile {outfile} exists, aborting. Use -t if you would like to use this table as input.\n".format(
                    outfile=outfile))
            sys.exit()

        # generate tab-delimited OTU matrix with taxonomy in final column
        # alterantive would be the Python biom functions
        command = f"biom convert -i {infile} -o {infile}.txt --to-tsv --header-key taxonomy"
        sys.stderr.write("Executing: {command}\n".format(command=command))
        os.system(command)

        infile = infile + '.txt'

    # read the tab-delimited data & collapse
    print(infile)
    sys.stderr.write("Reading input data\n")

    collapsed, taxonomy, samples = read_OTU_table(infile)
    taxon_to_trace = traces_from_taxonomy(collapsed)
    trace_to_taxon = inverse_dict(taxon_to_trace)
    counts = infer_internal_counts(collapsed)

    # printing the results
    sys.stderr.write("Printing output\n")
    traces = list(trace_to_taxon.keys())
    traces.sort()
    samples.sort()
    
    outtext = [("#class\tclass id\t" + "\t".join(samples))]

    for trace in traces:
        line = [trace]
        nice_taxon = trace_to_taxon[trace][-1]
        nice_taxon = nice_taxon.split('__')
        if nice_taxon[0] == 'Unknown':
            nice_taxon = "Unknown"
        else:
            nice_taxon = label_replace[nice_taxon[0]] + " - " + nice_taxon[-1]
        line.append(nice_taxon)
        taxon = trace_to_taxon[trace]
        for sample in samples:
            line.append("%f" % (counts[taxon][sample]))
        outtext.append("\t".join(line))

    if options['outfile'] == 'stdout':
        print('\n'.join(outtext))
    else:
        output = open(options['outfile'], 'w')
        output.write('\n'.join(outtext))
        output.close()
