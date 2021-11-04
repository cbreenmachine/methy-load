#!/bin/bash
input_dir="../../Fastq/" 

ls -1 ${input_dir}*.fastq.gz | sed -E s/[1-2]_001.fastq.gz// | uniq \
	| parallel -j 12 --gnu --env PATH 'trim_galore --fastqc --phred33 --output_dir ../../FastqTrimmedPair/ --paired {}1_001.fastq.gz {}2_001.fastq.gz'
