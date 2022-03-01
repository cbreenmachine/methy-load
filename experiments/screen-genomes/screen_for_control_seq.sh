#!/bin/bash
# Runs alignment on a few genomes.
# See fastq screen here: https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/
# Dont (cant) run pair ended: https://github.com/StevenWingett/FastQ-Screen/issues/32
input_dir="../../data/2021-11-03-batch01/pool01/group1/01-fastq-trimmed/"

ls -1 ${input_dir}100*fq | uniq > PREFIX

for P in `cat ./PREFIX`;
do
fastq_screen --conf ../../../FastQ_Screen_Genomes_Bisulfite/fastq_screen.conf \
    --tag --aligner bowtie2 \
    --bisulfite ${P}
done
    