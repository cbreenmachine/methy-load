
# wgbs-load

This repository contains the scripts used to align and otherwise process whole-genome bisulfite sequenced data collected on ~300 patients: one third with late-onset Alzheimer's disease, one-third with mild cognitive impairment, and one-third healthy.

# Organization

Data is coming in pools of five samples. Possible batch effects?


Trimmed files keep the structure within
```
wgbs-load
|--data
  |--2021-01-01-batch1
  |  |--pool01
  |  |  |--00-fastq
  |  |  |  |--sample01_R1.fq
  |  |  |  |--sample01_R2.fq
  |  |  |--01-fastq-trimmed
  |  |  |  |--sample01_R1-trimmed.fq
  |  |  |  |--sample01_R2-trimmed.fq
  |  |  |--02-mapped-reads
  |  |  |--03-called
  |  |--pool02
  |  |--00-fastq-trimmed
  |  |   |--sample01_R1_trimmed.fq
  |  |   |--sample02_R1_trimmed.fq
```


# Environment

As much as possible, 

# Methods

Six samples were sent to University of Illinois for whole-genome bisulfite sequencing. They used the EM kit (EpiMark?).

EpiMark.
