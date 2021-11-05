
# wgbs-load

This repository contains the scripts used to align and process whole-genome bisulfite sequenced data collected on ~300 patients: one third with late-onset Alzheimer's disease, one-third with mild cognitive impairment, and one-third healthy.

# Organization

Data is coming in "pools" of five samples. The structure of the data directory reflects how we received the data. As of November 4, 2021, we've recieved 20 samples (four pools of five).

This is also tracked in the


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

There are three "tiers" of data that we run throught the pipeline. The first, located in `./data/dummy/toy-data/` has two samples, each with forward and reverse reads. In reality, `dummy1_R[1-2].fq (dummy2_R[1-2].fq)` is just the first (last) 250,000 reads from sample 110. These files are used to debug the `code/runall.sh` pipeline. Trimming takes a few seconds, and mapping/calling are also quick. The second dataset, located in `./data/dummy/pool1/` more faithfully replicates the structure/naming of the actual data. A copy of sample 110's raw reads is located in this folder. To run this pipeline start to finish on one sample take about one day (24 hours), although running five or so samples should take aout the same amount of time due to parallelization. Currently, the bottleneck is trimming. `cutadapt` and `trim_galore` have 

# Environment

As much as possible, we use one conda environment, which can be created on your device with
`conda create --name load --file=load.yml`

Most code is run on a cluster of servers which run [Scientific Linux](https://en.wikipedia.org/wiki/Scientific_Linux) 7.9 (Nitrogen), which is based on Red Hat.

Environments are managed through [MiniConda](https://docs.conda.io/en/latest/miniconda.html).

# Methods

To process we use the [gemBS pipeline](https://github.com/heathsc/gemBS-rs), with some checking with [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). We also use [GNU parallel](https://www.gnu.org/software/parallel/) whenever possible.