#! /bin/bash
parallel --jobs 4 Rscript 00-impute-missing-vals.R --chr {} ::: chr6 chr18
parallel --jobs 4 Rscript 01-find-DMRs.R --DSS --chr {} ::: chr6 chr18
