#!/usr/bin/env python

# extract_methylation.py
# Takes .bcf output from gemBS pipeline and outputs a tsv

# Data
# Here's one line from an example bcf file (split over three lines):
# chrZZ    POS   .       G       .       QUAL      PASS    CX=TTGGT        
# GT:FT:DP:MQ:GQ:QD:GL:MC8:AMQ:CS:CG:CX   
# Some_other_stuff:     0,0,13,0,4,0,0,0    :some_other_stuff
# These eight numbers are 

# In MC8
# If the reference nt is 'C', then the number of unmethylated reads are #T (methylated is #C)
# If the reference nt is 'G', then the number of unmethylated reads are #A (methylated is #G)

# Considerations
# Coordinate system
#   1. vcf/bcf files use a 1-based coorinate system with inclusive bounds.
#       this means that a CpG on the first two positions would have POS
#   2. .BED files use 0-based 

# Output:
#   A (CSV or TSV or BED file)? Start with CSV and then move from there


#TODO include C context (currently only considering CpG on )
#TODO multiple output file supports

import argparse
import time
from pysam import VariantFile


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
    bcf_in = VariantFile("101.bcf")  # auto-detect input format
# The things in rec.format
# GT FT DP MQ GQ QD GL MC8 AMQ CS CG CX

for rec in bcf_in.fetch('chr1', 1000000, 1100000):
    d = get_field_dic(rec)
    if d['CG'] == 'Y':
        methylated, unmethylated = get_methylation_estimate(d['MC8'], d['CS'])
        print(methylated / (methylated + unmethylated))
        #TODO handle division by zero--just means there's a CpG site with no coverage!!
        #TODO coverage extracted from non-informative 
        #TODO make methylation function just accept dictionary as input...?
    
    
        