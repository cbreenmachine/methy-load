#!/bin/bash

idir="../../data/cov_meth/"
rm ../../data/prin_comps/var_explained.tsv

ls ${idir} > INPUT
parallel --link --workdir . --jobs 6 python compute_PCs.py --ifile ${idir}{1} :::: INPUT
rm INPUT