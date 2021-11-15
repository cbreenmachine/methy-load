#!/usr/bin/env python

# extract_methylation.py
# Takes .bcf output from gemBS pipeline and outputs a tsv

# In MC8
# If the reference nt is 'C', then the number of unmethylated reads are #T (methylated is #C)
# If the reference nt is 'G', then the number of unmethylated reads are #A (methylated is #G)

# Considerations
# Coordinate system
#   1. vcf/bcf files use a 1-based coorinate system with inclusive bounds.
#       this means that a CpG on the first two positions would have POS=(1,2)
#   2. BED files are zero-based so would have coordinates POS=(0,2)
# Output:
#   A (CSV or TSV or BED file)? Start with CSV and then move from there

#TODO coverage extracted from non-informative 
#TODO make methylation function just accept dictionary as input...?
#TODO coordinates to output file (may need zero or one based correction)
#TODO does DMRseq want coverage informative or coverage non-informative?
#TODO allow filtering by chromosome?  
#TODO include C context (currently only considering CpG on )
#TODO multiple output file supports (CSV, ENCODE style BED)

import argparse
import time
from pysam import VariantFile
from itertools import zip_longest
import csv


def get_field_dic(record):
    """returns the stuff we want in bcf line as a dictionary
    Input: one line (record) of vcf or bcf
    Output: dictionary with keys like 'MC8', 'GQ'
    """
    d = {}
    # quirky thing with pysam
    for sample in record.samples.items():
        for field in sample[1].items():
            # field is a tuple where the first item is the flag (e.g. MC8)
            # and the second item is the data we want: (0, 0, ....)
            d[field[0]] = field[1]
    return d

def get_methylation_estimate(mc8, strand):
    if strand == '+':
        # Number of Cs since Cs remained in tact if methylated
        methylated = mc8[5]
        # Number of Ts
        unmethylated = mc8[7]
    elif strand == '-':
        methylated = mc8[6] # Number Gs
        unmethylated = mc8[4] # Number As
    return(methylated, unmethylated)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Takes a list of input files? Or Idrectory...TBD")
    parser.add_argument("--input_file", default = "./101.bcf")
    parser.add_argument("--output_dir", default = "./extract_output/")
    parser.add_argument("--merge_strands", action = "store_true")

    args = parser.parse_args()

    infile = VariantFile("101.bcf", threads=4)
    csv_out_name = args.input_file.replace('.bcf', '.csv')
    ofile = open(csv_out_name, "w")

    # Column names for ouptut
    writer = csv.writer(ofile)
    writer.writerow(["chr", "pos", "reference", "call", "methylated", "unmethylated", "strand"])

    # The things in rec.format
    # GT FT DP MQ GQ QD GL MC8 AMQ CS CG CX
    # 480 minutes per one bcf--unacceptable!!!
    # 7 minutes for chrom 22--using 4 threads
    # 7 minutes for chrom 22--using 8 threads

    # Iterator
    #I = infile.fetch('chr1', 100000, 110000)
    I = infile.fetch('chr22')

    # Iterate two records at a time if merging...
    #for rec1, rec2 in zip_longest(*[I]*2):
    for rec2 in I:
        #data_1 = rec1.samples.items()[0][1].items()
        data_2 = rec2.samples.items()[0][1].items()

        # rec2 should be the base. Is it CpG? Then do the conditional tests
        # rec2 can be negative strand (and should still be written out)
        if (data_2[10][1] == 'Y'):
            m2, um2 = get_methylation_estimate(data_2[7][1], data_2[9][1]) 
            # This is the merge condition. The records need to be one position away from each other
            # They need to both be CpGs
            # They need for the first position on 
            #if (rec2.pos - rec1.pos == 1 and data_1[10][1] == 'Y' and data_1[9][1] == "+" and data_2[9][1] == "-"):
            #    m1, um1 = get_methylation_estimate(data_1[7][1], data_1[9][1]) 
                # Pack the information we want to write out into a tuple
                #TODO how to record genotype likelihood?
             #   cleaned = [rec1.contig, rec1.pos, rec1.ref, data_2[0][1], m1 + m2, um1 + um2, '+-']
            #else:
            cleaned = [rec2.contig, rec2.pos, rec2.ref, data_2[0][1], m2, um2, data_2[9][1]]

            # Open/create writer and then write to output file
            writer = csv.writer(ofile)
            writer.writerow(cleaned)
            last_cleaned = cleaned
    ofile.close()