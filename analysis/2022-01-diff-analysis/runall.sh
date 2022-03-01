#! /bin/bash
parallel --jobs 4 Rscript 00-impute-missing-vals.R --chr {} ::: chr{1..22}
# parallel --jobs 4 Rscript 01-find-DMRs.R --DSS --chr {1} --smoothing {2} --num_pcs {3} ::: chr6 chr18 ::: 100 150 200 ::: 0 2 4
parallel --jobs 6 Rscript 01-find-DMRs.R --DSS --chr {1} --smoothing {2} --num_pcs {3} ::: chr{1..22} ::: 150 ::: 2

# Defaults
# parallel --jobs 6 Rscript 01-find-DMRs.R --DSS --chr {} --smoothing 100 --num_pcs 2 ::: chr6 chr18


>CHR && >EXP
for v in $(find . -name models.RData | grep PCs-2/);
do
    echo ${v}
    Rscript 02-plot-array-comps.R --ifile "${v}"
    chr=$(echo ${v} | cut -d "/" -f2) # >> CHR
    exp=$(echo ${v} | cut -d "/" -f3) # >> EXP

    echo ${chr} >> CHR
    echo ${exp} >> EXP
    Rscript 02-make-track-plots.R --chr "${chr}" --experiment "${exp}"
done


# parallel Rscript 03-count-genes.R --idir {} ::: ./chr6 ./chr18