#!/bin/bash

path="../../data/2021-11-01-illinois-batch1/00-fastq-110-only/"


zcat ${path}110_TGTTCGCC-TCCTACCT_L00M_R1_001.fastq.gz | head -4000000 > ${path}dummy1_R1.fq
zcat ${path}110_TGTTCGCC-TCCTACCT_L00M_R2_001.fastq.gz | head -4000000 > ${path}dummy1_R1.fq

zcat ${path}110_TGTTCGCC-TCCTACCT_L00M_R1_001.fastq.gz | tail -4000000 > ${path}dummy2_R1.fq
zcat ${path}110_TGTTCGCC-TCCTACCT_L00M_R2_001.fastq.gz | tail -4000000 > ${path}dummy2_R1.fq