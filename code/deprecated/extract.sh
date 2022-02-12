#! /bin/bash
ifile="$1"
ofile="$2"

bcftools query --format '%CHROM\t%POS\t%REF\t[%CG]\t[%CS]\t[%MC8]\n' --include 'DP>2 && CG="Y" && DP < 50' ${ifile} \
    | python extract.py --ofile ${ofile} --merge