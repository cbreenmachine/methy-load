#!/bin/bash
# SPecifically for 
indir="../../FastqTrimmedPair/"
out="../../meta.csv"
echo "Barcode,Dataset,End1,End2" > $out

for f in $indir*1.fq.gz;
do
    barcode="s"$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/RA//')
    dataset=$(echo ${f##*/} | sed -E 's/_[ACTG-]{1,20}_L002//' | sed -E 's/val_[1-2].fq.gz//')
    file1=$(echo ${f##*/})
    file2=$(echo ${f##*/} | sed -E 's/R1_001/R2_001/')
    echo $barcode","$dataset","$file1","$file2 >> $out
done
cat $out