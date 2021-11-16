#!/usr/bin/env python
import csv
import sys
import itertools


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

x = 0
y = 0
z = 0

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
            x += 1
        elif skip_me:
            skip_me = False
            y += 1
        else:
            my_writer.writerow(fla)
            z += 1
    except:
        my_writer.writerow(a.split() )

print(x,y,z)