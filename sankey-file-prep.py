#!/usr/bin/python3

# Author: Harm Laurense
# Last edited: 15-08-2025 [by TE]
# Function: This script is used to create the necessary .csv files to create sankey diagram(s) using R.
#
# Parameter 1: Integer/Float for filtering of taxa based on their (average) rel. abund.
# Parameter 2: String ( boolean principle* ) for creating a .csv file for each individual sample.
# Parameter 3: String input for the metadata file.
# Parameter 4: String ( booleon principle* ) for creating a .csv file for each unique rankstat combination of samples.
# Parameter 5: String input for the biotaviz file.
# * checks if input is equal to "true"
#
# Typical run:
# python3 .\sankey-file-prep.py 0.01 false false metadata.tsv relative-table.biotaviz.txt

import csv
import sys
import traceback
import itertools
import argparse

# Dictionary used to determine the numbers for linking nodes (based on taxonomic rank)
taxonomic_ranks_dict = {
    "empty": -1,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
    "species": 6}
# Dictionary used to determine the numbers for linking nodes (based on taxonomic rank / last used taxonomic rank)
taxonomic_ranks_last_linked_rank_dict = {
    "phylum": 0,
    "class": 0,
    "order": 0,
    "family": 0,
    "genus": 0,
    "species": 0}

def main(tax_filter, sample_repeat, mappingfile, combine_rankstat):
    """
    Determines what functions are needed to be called based on command line input.
    :param mappingfile: File containing the metadata. Determines which samples are averaged.
    :param tax_filter: Parameter for filtering (low) relative abundance.
    :param sample_repeat: Parameter which determines if files are created for every individual sample.
    :return: .csv files according to user input, to be used in R script for creating the sankey diagrams
    """
    # [DEFAULT] Create sample average file over all samples (includes samples without metadata values)
    sample_average_all()

    # [REPEAT = TRUE] Create a .csv file for every sample (needed for generating sankey diagram in R script)
    if sample_repeat.lower() == "true":
        total_samples = determine_sample_total()
        for sample in range(total_samples):
            hierarchy_counts(tax_filter, sample, average_all_samples, filename_combination)

    # [DEFAULT] Create sample average file each Rankstat column
    all_sets = get_sets(mappingfile)
    if len( all_sets.keys() ) > 0 :
		
        for set in all_sets :
            filename_rankstatheaders, rankstat_samples = determine_rankstat_samples(set, all_sets[set])
            sample_index = determine_sample_index()				
            indexed_combinations = combination_to_index(sample_index, rankstat_samples)

            allcolumnsamples = sum(indexed_combinations, [])
            sample_average(allcolumnsamples, set) # for all samples, in one Rankstat column

            # [COMBINE = TRUE] Create sample average file for every individuel study group in a Rankstat column
            if combine_rankstat.lower() == "true": 
                for index, combination in enumerate(indexed_combinations):
                    sample_average(combination, filename_rankstatheaders[index]) # for samples in each different study group, in one Rankstat column

def determine_sample_total():
    """
    Determine the total amount of samples by counting the columns. The first 2 columns aren't samples and thus skipped.
    :return:Number of samples (total)
    """
    try:
        with open(biotavizfile, "r") as file:
            line = file.readline()
            total_samples = len(line.rstrip().split('\t')[2:])
        return total_samples
    except IndexError:
        print("# IndexError; check if the correct file is given as input: ", traceback.print_exc())

def hierarchy_counts(tax_filter, sample, average_samples, filename_combination):
    """
    Generate the files necessary for creating a sankey diagram from the biotaviz file.
    :param tax_filter: Parameter for filtering (low) relative abundance.
    :param sample: Integer (standard 0) used as index. This only changes if sample_repeat is set to true.
    :param average_samples: List of average values (currently from all rankstat sample combinations).
    :param filename_combination: Specific string correlating to the unique combination, needed for unique filenames.
    :return: Variables (link1, link2, label) are determined and finally given to the write_new_biotaviz_file() function.
    Which in return will write the necessary .csv files.
    """
    link1 = ["link1", "link1"]
    link2 = ["link2", "link2"]
    label = [["label", "value"]]
    count = 0
    nonzeros = 0
    last_rank = "empty"
    removed_entries = []
 
    try:
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            taxonomic_rank = []
            for index, line in enumerate(file.readlines()):
                if line.strip():
                    line = line.rstrip().split('\t')
                    if sample == 'AVRG':
                        tax_value = float(average_samples[index])
                        if tax_value > 0 : nonzeros += 1
                    else:
                        # Skip first 2 columns
                        tax_value = float(line[sample + 2])
                        if tax_value > 0 : nonzeros += 1
                   
                    # Values of 0 (relative abundance) or below taxonomic filter (standard 1%) aren't used
                    if tax_value >= tax_filter and tax_value > 0:
                        tax_rank = line[1].split('-')[0].rstrip()
                        tax_specific = line[1].split('-')[1].rstrip()
                        taxonomic_rank.append(tax_rank)
                        tax_with_score = tax_specific + ":" + str(round(float(tax_value) * 100, 2)) + "%"
                        # This var replacement was used during testing for various value sizes based on taxonomic rank
                        # tax_value = taxonomic_ranks_valuesize_dict.get(tax_rank)*tax_value
                        label.append([tax_with_score, tax_value])
                    else:
                        removed_entries.append(line)

            # The following is the logic to determine the number combinations for connecting the nodes
            for rank in taxonomic_rank:
                if rank == "domain":
                    continue
                elif taxonomic_ranks_dict.get(rank):
                    count += 1
                    if rank == last_rank:
                        link1.append(taxonomic_ranks_last_linked_rank_dict.get(rank))
                    if taxonomic_ranks_dict.get(last_rank) < taxonomic_ranks_dict.get(rank):
                        if link2[-1] != "link2":
                            link1.append(link2[-1])
                        else:
                            link1.append(0)
                    elif taxonomic_ranks_dict.get(last_rank) > taxonomic_ranks_dict.get(rank):
                        link1.append(taxonomic_ranks_last_linked_rank_dict.get(rank))
                    link2.append(count)
                    last_rank = rank
                    taxonomic_ranks_last_linked_rank_dict.update({rank: link1[-1]})

        if nonzeros == 0 :
            print("# WARNING: The following sample contains only zero values, we will therefore skip the creation of a Sankey plot for :", sample) 
        elif len(label) > 1:
            write_new_biotaviz_file(link1, link2, label, sample, average_samples, filename_combination)
        else:
            sys.exit(print("# No matches with current criteria found, try lowering the given filter for relative abundance"))
    except IndexError:
        print("# IndexError; check if the correct file is given as input: ", traceback.print_exc())

def write_new_biotaviz_file(link1, link2, label, sample, average_samples, filename_combination):
    """
    Write a .csv file containing the necessary information for creating a sankey diagram (used in: Sankey R module)
    :param link1: List of numbers which represents the node being connected from.
    :param link2: List of (second) numbers which represents the node connected to.
    :param label: List of lists containing the labels (taxonomic rank : % abundance) and value (relative abundance)
    :param sample: Integer (index) used to create an unique filename
    :param average_samples: Needed to determine which filename is to be used.
    :param filename_combination: Specific string correlating to the unique combination, needed for unique filenames.
    :return: .csv file (4 columns)
    """
    try:
        if not average_samples:
            filename = f"biotaviz_sankey_prepfile-{sample}.csv"
        else:
            filename = f"biotaviz_sankey_prepfile-{filename_combination}.csv"
        with open(filename, "w", newline="") as f:
            wr = csv.writer(f)
            # column headers
            for index, number in enumerate(link1):
                wr.writerow([number, link2[index], label[index][0], label[index][1]])
    except IndexError:
        print("# IndexError; can't write to biotaviz_sankey_prepfile-{sample}.csv: ", traceback.print_exc())

def get_sets(filename):
    """
    Determine what subsets should be compared (rankstat)
    :param filename: Metadata file (should contain rankstat columns)
    :return: Dictionary containing each rankstat column, containing all unique values with corresponding samples
    """
    lines = load_txt(filename).strip().split('\n')
    headers = {}
    datasets = {}
    for i, element in enumerate(lines[0].split('\t')):
        if element.split('_')[0].lower() == 'rankstat':
            dataset = element.split("_")[1]
            headers[i] = dataset
            datasets[dataset] = {}
    for line in lines[1:]:
        lineg = line.split('\t')
        sample_id = lineg[0]
        for i, element in enumerate(lineg):
            if i in list(headers.keys()) and element != '':
                sample_class = element
                dataset = headers[i]
                if sample_class not in datasets[dataset]:
                    datasets[dataset][sample_class] = []
                datasets[dataset][sample_class].append(sample_id)
    return datasets

def load_txt(filename):
    """
    Opens files and returns their input/content
    :param filename: File
    :return: File content
    """
    fileinput = open(filename)
    text = fileinput.read()
    fileinput.close()
    return text

def determine_sample_index():
    """
    Determines the index for each rankstat column, necessary in following functions to determine the sample average
    :return: Dictionary containing the rankstat columns and their corresponding index.
    """
    try:
        sample_index = {}
        with open(biotavizfile, "r") as file:
            sample_headers = file.readline()  # skip column headers + root
            sample_headers = sample_headers.rstrip().split("\t")
            for index, header in enumerate(sample_headers[2:]):
                sample_index.update({header: index + 2})
        return sample_index
    except IndexError:
        print("# IndexError; check if the correct file is given as input: ", traceback.print_exc())

def determine_rankstat_samples(header, set):
    """
    Create two lists necessary for unique filenames and to determine index + average of samples in following functions
    :param all_sets: Dictionary containing each rankstat column, containing all unique values with corresponding samples
    :return: A list containing the filenames (rankstat+unique_value) and a list with their corresponding samples
    """
    rankstat_samples = []
    filename_rankstatheaders = []
    for key in set:
        column_value = header + "-" + key
        filename_rankstatheaders.append(column_value)
        sample_value = set[key]
        rankstat_samples.append(sample_value)

    return filename_rankstatheaders, rankstat_samples

def combination_to_index(sample_index, unique_sample_combinations):
    """
    Create a nested list with the index for each sample of each combination needed for determining the avarage in
    following functions
    :param sample_index: Dictionary containing the rankstat columns and their corresponding index.
    :param unique_sample_combinations: Nested list containing their corresponding (unique/no duplicate) samples
    :return: Nested list containing the index for each sample (for each unique combination)
    """
    indexed_combinations = []
    for index, combination in enumerate(unique_sample_combinations):
        indexed_combinations.append([])
        for unique_sample in combination:
            indexed_combinations[index].append(sample_index.get(unique_sample))
        indexed_combinations[index] = sorted(indexed_combinations[index])
    return indexed_combinations

def sample_average(combination, combination_header):
    """
    Calculates the average values of samples for each combination of rankstat column values.
    :param combination: Unique combination of rankstat column values
    :param combination_header: Unique filename for each combination
    :return: .csv file (used to create sankey for sample average)
    """
    try:
        average_samples = []
        rankstat_values = []
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            for line in file:
                if line.strip():
                    line = line.rstrip().split('\t')
                    for rankstat_column_index in combination:
                        rankstat_values.append(float(line[rankstat_column_index]))
                    average_samples.append(sum(rankstat_values) / len(rankstat_values))
                    rankstat_values = []
        if ( all(i == 0 for i in average_samples) ) :
            print("# WARNING: The following combination of study groups contains only zero values, we will therefore skip the creation of a Sankey plot for :", combination_header)
        else :
            hierarchy_counts(tax_filter, 'AVRG', average_samples, combination_header)

    except IndexError:
        print("# IndexError; check if the correct file is given as input: ", traceback.print_exc())

def sample_average_all():
    """
    Calculates the average values of all samples, to be used in generating a .csv file. (even without value in rankstat)
    :return: .csv file (used to create sankey for sample average)
    """
    try:
        average_samples_all = []
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            for line in file:
                if line.strip():
                    line = line.rstrip().split('\t')
                    all_samples = []
                    for value in line[2:]:
                        all_samples.append(float(value))
                    average_samples_all.append(sum(all_samples) / len(all_samples))
        filename_part = "AverageAllSamples"

        if ( all(average_samples_all) == 0 ) :
            sys.exit(print("# ERROR: The average of all samples contains only zero values, we will therefore skip the creation of a Sankey plot for : AverageAllSamples"))
        else : 
            hierarchy_counts(tax_filter, 'AVRG', average_samples_all, filename_part)

    except IndexError:
        print("# IndexError; check if the correct file is given as input: ", traceback.print_exc())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sankey file preparation", add_help=True)
    parser.add_argument('--taxa-filter', dest='tax_filter', help='Taxa filter', default=0.01, type=float)
    parser.add_argument('--sample-repeat', dest='sample_repeat', help='Sample repeat, default is false', default="false")
    parser.add_argument('--combine-rankstat', dest='combine_rankstat', help='Combine rankstat, default is false', default="false")
    parser.add_argument('-i', dest='infile', help='name of input file', required=True)
    parser.add_argument('-m', dest='mapping', help='name of mapping file', required=True)
    options = vars(parser.parse_args())

    # Global variables
    filename_combination = ""
    biotavizfile = options['infile']
    ## sample = "AVRG"
    average_all_samples = ""
    tax_filter = options['tax_filter']
    sample_repeat = options['sample_repeat']
    mappingfile = options['mapping']
    combine_rankstat = options['combine_rankstat']

    # Main
    if not (tax_filter > 0 and tax_filter < 1):
        sys.exit(print("# Use a number between 0 and 1 as parameter for filtering relative abundance"))
    try:
        main(tax_filter, sample_repeat, mappingfile, combine_rankstat)
    except ValueError:
        print("# Parameter given was not a valid numeric value: ", traceback.print_exc())
        print("# If the input is a decimal number, use a decimal point instead of comma (eg 0.01 instead of 0,01)")
    except IndexError:
        print("# Not enough parameters were given: ", traceback.print_exc())
