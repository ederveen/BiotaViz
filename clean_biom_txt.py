#!/usr/bin/env python3
import sys
from argparse import ArgumentParser

# settings
undefined_labels_part = ["unclassified", "uncultured", "unknown", "Unclassified", "Uncultured", "unknown", "metagenome"]
undefined_labels_full = ["_", "", "__"]
description = "Add last-known level to biom taxon trace"


def load_txt(infile):
    file_input = open(infile, 'r')
    lines = file_input.read().rstrip().split('\n')
    file_input.close()
    return lines

def clean_trace(tax_trace):
    last_known = tax_trace[0]
    for i, taxon in enumerate(tax_trace):
        undefined = False
        empty = False
        for element in undefined_labels_part:
            if taxon.find(element) != -1:
                undefined = True
        for element in undefined_labels_full:
            clean_taxon = taxon[3:] 
            if clean_taxon == element:
                empty = True
                undefined = True
        if undefined is True:
            if empty is True:
                taxon += "Unclassified"
            # new_taxon = taxon + '_'+last_known.replace("__", "_")
            # tax_trace[i] = new_taxon
            tax_trace[i] = ''
        else:
            last_known = taxon
    return tax_trace


# main program
if __name__ == '__main__':

    parser = ArgumentParser(description=description, add_help=True)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-o', dest='outfile', help='name of output file', required=True)
    options = vars(parser.parse_args())

    infile = options['infile']
    outfile = options['outfile']

    lines = load_txt(infile)
    if lines[1].split('\t')[-1].lower() != "taxonomy":
        sys.stderr.write('Last header of second line is not "taxonomy", aborting\n')
        sys.exit()
    else:
        taxonomy_column = len(lines[1].split('\t'))-1

    output = open(outfile, "w")
    output.write(lines[0]+'\n'+lines[1]+'\n')
    for line in lines[2:]:
        lineg = line.split('\t')
        tax_trace = lineg[taxonomy_column].split('; ')
        tax_trace = clean_trace(tax_trace)
        tax_trace = '; '.join(tax_trace)
        new_line = '\t'.join(lineg[:taxonomy_column] + [tax_trace])
        output.write(new_line+'\n')
    output.close()
