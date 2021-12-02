#!/bin/bash
# Just import the functions...
import csv
import sys
import itertools
import argparse

#############################################################
###### Get (unmerged) methylation numbers for CpG sites #####
#############################################################
def get_methylation_estimate(mc8, strand):
    '''manipulates 
    '''
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
    # 0:3 are chromosome, position, ref, CpG (keeping for when we have CHH, etc)
    # 4 becomes +- to indicate it's merged
    # TODO coverage
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



#############################################################
############# Merge CpG sites and record ####################
#############################################################


# https://stackoverflow.com/questions/5434891/iterate-a-list-as-pair-current-next-in-python
def pairwise(iterable):
    """generates two iterators, second offset by one
    Input: file iterator
    Output: prints first line of iterator
    """
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)   
    return zip(a, b)   

def read_my_line(a):
    tmp = a.split() # coming in as tab-separated
    tmp[1] = int(tmp[1])
    tmp[5] = int(tmp[5])
    tmp[6] = int(tmp[6])
    return tmp


# Just to make it more portable/flexible if we change I/O
infile = sys.stdin
outfile = sys.stdout

# Keep tab delimited, easier to read and takes the same amount of storage
my_writer = csv.writer(outfile, delimiter='\t')
skip_me = False

for a, b in pairwise(infile):

    try:
        fla = read_my_line(a) # field list a
        flb = read_my_line(b)

        # Possibilty that we have consecutive positions
        if (flb[1] - fla[1] == 1 and fla[4] == "+"):
            # Take all the usual information, but add methylated reads and unmethylated reads together
            # and call the strand +- to mean "merged"
            my_writer.writerow( fla[0:4] + ["+-", fla[5] + flb[5], fla[6] + flb[6] ])
            skip_me = True 
        elif skip_me:
            skip_me = False
        else:
            my_writer.writerow(fla)
    except:
        my_writer.writerow(a.split() )



#############################################################
############# Loops through ...##################
#############################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='methylation extraction options')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
    parser.add_argument('--crawl', dest='accumulate', action='store_true',
                    help='should the script crawl (recursively)')

args = parser.parse_args()
    # want argparse--incorporate crawler
    # accept input output
    # accept --merge
    # accept --crawl (print files that will be input and output--store as variables not text)
    # 