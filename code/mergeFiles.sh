#!/bin/bash
idir="../../FastqTrimmed/"

# RA110_CACTGTAG-AGTCGCTT_L002_R1_001.fastq.gz
# RA110_CACTGTAG-AGTCGCTT_L002_R2_001.fastq.gz

ls -1 ${idir}*.fq.gz | rev | cut -d "/" -f1 | rev | sed -E s/_L002_R[1-2]_001_trimmed.fq.gz// | uniq \
    | parallel -j 6 --gnu --env PATH "cat "$idir"{}_L002_R1_001_trimmed.fq.gz "$idir"{}_L002_R2_001_trimmed.fq.gz > ../../FastqMerged/{}_trimmed_merged.fq.gz"