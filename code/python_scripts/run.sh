#!/bin/bash

idir="../../data/cov_meth/"

ls ${idir} > INPUT
parallel --link --workdir . --jobs 6 python compute_PCs.py --ifile ${idir}{1} :::: INPUT
rm INPUT