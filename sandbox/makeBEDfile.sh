#!/bin/bash

# bcftools view ../../PilotDHMRI/calls/sample006_112.bcf | head -300 | makeBEDfile.sh -o test.bed

# makeBEDfile.sh extracts the most in
# takes two command line arguments -i (input file should be .vxf or .bcf) and -o should be .csv
# Runs like `makeBEDfile.sh -i ../../PilotDHMRI/calls/ -o ../../DataRaw/BED_CpGs`

# Handle -i and -o flags
#while getopts i:o: flag
#do
#    case "${flag}" in
#        i) idir=${OPTARG};;
#        o) odir=${OPTARG};;
#    esac
#done

#mkdir -p "${odir}"

### Workhorse function
makeBED () {

out="$(echo $1 | sed -E s/bcf/bed/)"
echo -e "chrom\tpos\trefBP\tref5\tcalled5\tCpG\tstrand\treaddepth\tA\tC\tG\tT\tfilterInfo\tmethylation\tnumerator\tdenominator" > $out
# Square brackets are used to extract FORMAT tags   
bcftools query \
	--format '%CHROM\t%POS\t%REF\t%CX\t[%CX]\t[%CG]\t[%CS]\t[%DP]\t[%MC8]\n' \
	--include 'DP>2 & CG="Y"' \
	--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
	$1 | awk -v OFS='\t' -f inferMethylation.awk >> $out 

}

makeBED 101.bcf
