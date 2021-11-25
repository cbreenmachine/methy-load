#!/usr/bin/env python
import csv
import sys

def get_methylation_estimate(mc8, strand):
    methylated = float('NaN')
    unmethylated = float('NaN')

    if strand == '+':
        # Number of Cs since Cs remained in tact if methylated
        methylated = mc8[5]
        # Number of Ts
        unmethylated = mc8[7]
    elif strand == '-':
        methylated = mc8[6] # Number Gs
        unmethylated = mc8[4] # Number As
    return(methylated, unmethylated)


def merge_two_lines(line1, line2):
    return(line1[0:3] + ["+-", line1[5] + line2[5], line1[6] + line2[6]])

infile = sys.stdin
outfile = sys.stdout

my_writer = csv.writer(outfile,delimiter='\t')
my_writer.writerow(["chr", "pos", "reference", "context", "strand", "methylated", "unmethylated"])

for line in infile:
    fields_list = line.split() # coming in as tab-separated
    mc8 = [int(x) for x in fields_list[5].split(",")] # grab read counts as list of integers
    m, um = get_methylation_estimate(mc8, fields_list[4])
    my_writer.writerow(fields_list[0:3] + ["CpG", fields_list[4]] + [m, um])
