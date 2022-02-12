#!/usr/bin/env python
import csv
import sys
import itertools
import argparse

# HELPER FUNCTIONS
def get_methylation_estimate(mc8, strand):
    '''given the mc8 field from vcf/bcf, tell me the number of methylated
    and unmethylated counts
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
    '''take two two lines from bcf file and combine them (merge strand)
    '''
    return(line1[0:3] + ["+-", line1[5] + line2[5], line1[6] + line2[6]])

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


# WORKHORSE FUNCTION TO EXTRACT FROM BCF
def extract(infile, outfile):
    '''workhorse function that takes input from bcf and
    outputs 
    '''
    print("Opening writer for extraction...")

    my_writer = csv.writer(outfile, delimiter=',')
    # Assumes certain columns from standard in
    #TODO: type checks here
    my_writer.writerow(["chr", "pos", "reference", "context", "strand", "methylated", "unmethylated"])

    for line in infile:
        fields_list = line.split() # coming in as tab-separated
        mc8 = [int(x) for x in fields_list[5].split(",")] # grab read counts as list of integers
        m, um = get_methylation_estimate(mc8, fields_list[4])
        my_writer.writerow(fields_list[0:3] + ["CpG", fields_list[4]] + [m, um])
    print("Finished formatting... Will merge if requested...")

def merge_strands(infile, outfile):

    my_writer = csv.writer(outfile, delimiter=',')
    skip_me = False 
    x, y, z = 0, 0, 0

    for a, b in pairwise(infile):

        try:
            fla = read_my_line(a) # field list a
            flb = read_my_line(b)

            # Handle the possibilty that we have consecutive positions
            if (flb[1] - fla[1] == 1 and fla[4] == "+"):
                # Take all the usual information, but add methylated reads and unmethylated reads together
                # and name the strand +- to mean "merged". May change to '*'
                my_writer.writerow( fla[0:4] + ["+-", fla[5] + flb[5], fla[6] + flb[6] ])
                skip_me = True 
                x += 1
            # Handle the possibility that 
            elif skip_me:
                skip_me = False
                y += 1
            else:
                my_writer.writerow(fla)
                z += 1
        except:
            my_writer.writerow(a.split() )


def main(infile, outfile, merge=False):
    if merge:
        print("Will merge strands")
        extract(infile, sys.stdout) # bridge the two by writing / reading to standard ()
        merge_strands(sys.stdin, outfile)
    else:
        extract(infile, outfile)


if __name__== "__main__":
    # argparsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--ofile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--merge', action='store_true', default=False)
    args = parser.parse_args()

    main(args.ifile, args.ofile, args.merge)