#! /bin/bash
#parallel --jobs 4 Rscript 00-impute-missing-vals.R --chr {} ::: chr6 chr18 
parallel --jobs 4 Rscript 01-find-DMRs.R --DSS --chr {1} --smoothing {2} --num_pcs {3} ::: chr18 ::: 0 250 500 ::: 0 10 20


# Defaults
Rscript 01-find-DMRs.R --DSS --chr chr18 --smoothing 100 --num_pcs 2
Rscript 01-find-DMRs.R --DSS --chr chr6 --smoothing 100 --num_pcs 2


>CHR && >EXP
for v in $(find . -name models.RData);
do
    Rscript 02-plot-array-comps.R --ifile "${v}"
    echo ${v} | cut -d "/" -f2 >> CHR
    echo ${v} | cut -d "/" -f3 >> EXP
done
parallel --link --jobs 4 Rscript 02-make-plots.R --chr {1} --experiment {2} :::: CHR :::: EXP
rm CHR EXP