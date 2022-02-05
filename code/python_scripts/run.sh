#!/bin/bash
idir="../../data/cov-meth/"

ls ${idir} > INPUT
parallel --link --workdir . --jobs 8 python compute_PCs.py --ifile ${idir}{1} --filter_samples :::: INPUT
rm INPUT
