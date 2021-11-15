#!/usr/bin/env python
import csv
import sys
import itertools

# https://stackoverflow.com/questions/5434891/iterate-a-list-as-pair-current-next-in-python
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
   
    return zip(a, b)   


my_writer = csv.writer(sys.stdout)
my_writer.writerow(["chr", "pos", "something"])

infile = sys.stdin
outfile = sys.stdout

my_writer = csv.writer(outfile, delimiter='\t')

skip_me = False

for a, b in pairwise(infile):
    fl1 = a.split() # coming in as tab-separated
    fl2 = b.split() # coming in as tab-separated

    try:
        # Possibilty that we have consecutive positions
        if not skip_me:
            if fl2[1] - fl1[1] == 1:
                my_writer.writerow(fl1[0:4] + ["+-", fl1[6] + fl2[6], fl1[7] + fl2[7]] )
                skip_me = True #
            else:
                my_writer.writerow(fl1)
    except: 
        my_writer.writerow(fl1)