#!/usr/bin/env python3
# todoc

"""
Biotaviz_counts_to_abundance
---------------
.. module:: Biotaviz_counts_to_abundance
  :synopsis: Convert raw-count BiotViz-style file to relative abundance version, new edition
.. moduleauthor:: Jos Boekhorst

Convert raw-count BiotViz-style file to relative abundance version, new edition

Typical run::

    Biotaviz_counts_to_abundance.py -i BiotaViz.txt -o BiotaViz_relative.txt

Run the script with '-h' for a list of options.
"""
from argparse import ArgumentParser
import sys

############
# SETTINGS #
############

def checksZeroDivision(num1, num2):
    if num1 == 0.0 or num2 == 0.0:
        return 0.0
    else:
        return num1 / num2

# note: default parameters are set in argparse object (top of __main__)
description = "Convert raw-count BiotViz-style file to relative abundance version. First non-header line must be the root (i.e., total count)!"

# main program
if __name__ == '__main__':
    parser = ArgumentParser(description=description, add_help=True)
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-o', dest='outfile', help='name of output file', default='', required=False)
    parser.add_argument('-r', dest='root_name', help='taxon to take as root (i.e., set to 1)', default='', required=False)
    options = vars(parser.parse_args())

    if options['outfile'] == "":
        options['outfile'] = options['infile'].replace('.txt', '_relative.txt')

    with open(options['infile'], 'r') as f:
        lines = f.read().rstrip().split('\n')

    with open(options['outfile'], 'w') as f:
        f.write(lines[0] + '\n')
        if options['root_name'] != "":
            found = 0
            for line in lines:
                lineg = line.split('\t')
                if lineg[1] == options['root_name']:
                    root_trace = lineg[0]
                    total_counts = [float(element) for element in lineg[2:]]
                    found = 1
                    break
            if found == 0:
                sys.stderr.write('Could not find specified root name "' + options['root_name'] + '"\n')
                sys.exit(1)
        else:
            total_counts = [float(element) for element in lines[1].split('\t')[2:]]
        for line in lines[1:]:
            lineg = line.split('\t')
            if options['root_name'] != "":
                if lineg[0][:len(root_trace)] != root_trace:
                    skip = 1
                else:
                    skip = 0
            else:
                skip = 0
            if skip == 0:
                new_line = [lineg[0], lineg[1]]
                for i, count in enumerate([float(element) for element in lineg[2:]]):
                    new_line.append(f"{checksZeroDivision(count, total_counts[i])}")
                f.write('\t'.join(new_line) + '\n')
    sys.exit()
