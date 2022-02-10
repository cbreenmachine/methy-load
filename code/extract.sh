#! /bin/bash
input="$1"
output="$2"

# Consider DP < 100
bcftools query --format '%CHROM\t%POS\t%REF\t[%CG]\t[%CS]\t[%MC8]\n' --include 'DP>2 & CG="Y"' ${input} \
    | python extract_methylation.py | python merge_strands.py  > ${output}