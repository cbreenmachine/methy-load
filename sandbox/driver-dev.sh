
ALL="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

bcftools query --format '%CHROM\t%POS\t%REF\t[%CG]\t[%CS]\t[%MC8]\n' --include 'DP>2 & CG="Y"' 101.bcf 
	| python extract.py | python merge_strands.py  > 101_merge_chr22.tsv