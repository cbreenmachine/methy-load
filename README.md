
# wgbs-load

This repository contains the scripts used to align and process whole-genome bisulfite sequenced data collected on ~300 patients: one third with late-onset Alzheimer's disease, one-third with mild cognitive impairment, and one-third healthy.

# Organization

Link incoming...

# Environment

As much as possible, we use one conda environment, which can be created on your device with
`conda create --name load --file=load.yml`

Most code is run on a cluster of servers which run [Scientific Linux](https://en.wikipedia.org/wiki/Scientific_Linux) 7.9 (Nitrogen), which is based on Red Hat.

Environments are managed through [MiniConda](https://docs.conda.io/en/latest/miniconda.html).

# Methods

To process we use the [gemBS pipeline](https://github.com/heathsc/gemBS-rs), with some checking with [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). We also use [GNU parallel](https://www.gnu.org/software/parallel/) whenever possible.

# References

Organization of this repository is based in part on (that bioinformatics paper).
