#!/bin/bash
# runall.sh
# conda environment: runs in wgbs (preferred for parallel cutadapt) or load
# Runs entire bioinformatics pipeline starting with paired-end data and ending with bcf files.
# Made parallel with GNU parallel

# Steps in pipeline
# 1. Trimming using _trim_galore and cutadapt
# 2. gemBS mapping
# 3. gemBS calling

# Practical considerations
# 1. gemBS is a high memory pipeline--stick to the higher memory servers
# 2. Benchhmarking to come.
# 3. Run on ../data/dummy/toy[1-2]/ first to make sure messages and logs are 
# behaving how they're supposed to. (structure is analogous to real data with 'pool' swapped for 'toy')

#TODO: Make default environment wgbs (after finishes running)
# CONSTANTS
RUN_SERVERS="nebula-5,nebula-6,nebula-7"
DATE_STR=$(date +"%y-%m-%d")
MY_HOME=$(pwd)

# SETUP PATHS
#TODO: Make this accept from getopts
#ROOT_DIR="group0{1..4}/"
POOL_SUB_DIRS="$(echo ../data/2021-11-03-batch01/pool05-group0{1..4}/)" # be sure to have trailing slash

FASTQ_PATH="00-fastq/"
FASTQ_TRIMMED_PATH="01-fastq-trimmed/"
CALLS_PATH="03-calls/"
EXTRACT_PATH="04-extract/"
META_OUT="./meta.csv"
CONF_OUT="./conf.conf" 

for pool in ${POOL_SUB_DIRS}; do
# Functions are run in batch/pool/group/
# This directory contains 00-fastq
  cd ${pool}
  echo "Found these input files in ${ROOT_DIR}${pool} :"
  ls 00-fastq
  echo "Creating ${FASTQ_TRIMMED_PATH} directory"
  mkdir -p -v "${FASTQ_TRIMMED_PATH}"
  printf "\n" 
  cd ${MY_HOME}
done

echo "To continue running, hit y"
read -n 1 k <&1
if [[ $k = n ]] ; then
  echo "Re-configure directories or delete files as needed!"
elif [[ $k = y ]] ; then
  echo "continuing"
fi


for pool in $POOL_SUB_DIRS;
do
cd ${pool}
echo "Working on ${pool}"
#meta_out="./$(echo ${pool%/}).meta.csv" # name of meta file used in gemBS
#conf_out="./$(echo ${pool%/}).conf" # name of congiguration file used by gemBS
#fastq_path=${ROOT_DIR}${pool}"00-fastq/"
#fastq_trimmed_path=${ROOT_DIR}${pool}"01-fastq-trimmed/"

# trim files (conditional on them not being trimmed yet)
# Consider piping file names to TMP and then using parallel afterwords
ls -1 ${FASTQ_PATH}*R1*.fastq.gz | uniq > LEFT
ls -1 ${FASTQ_PATH}*R2*.fastq.gz | uniq > RIGHT

# 1. After trimming, rename files, then move them to appropriate dir
# 2. Make conditional ()
if [ -z "$(ls -A ${FASTQ_TRIMMED_PATH})" ]
then
   echo "${FASTQ_TRIMMED_PATH} is empty, trimming files now."
    # --link creates a mapping between the lines in LEFT and lines in RIGHT 
    # (one-to-one instead of pairwise combinations)
    # the fourth ':' means cat LEFT and RIGHT (don't treat as variable/expansion)
    parallel --link -S ${RUN_SERVERS} --workdir . --joblog ${DATE_STR}-trim.log \
        trim_galore --phred33 --cores 6 --output_dir ${FASTQ_TRIMMED_PATH} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT
else
   echo "${FASTQ_TRIMMED_PATH} is full (no need to trim/fastqc files); delete if you need to re-trim!"
fi

rm LEFT RIGHT
# END TRIMMING SECTION


# BEGIN MAKE META FILE
echo "barcode,dataset,end1,end2" > ${META_OUT}
for f in "${FASTQ_TRIMMED_PATH}"*val_1.fq
do
    # f="../../data/something.csv"
    # --> ${f##*/} is something.csv
    barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
    dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq//' | sed -E 's/_R1_//')
    file1=$(echo ${f##*/})
    file2=$(echo ${f##*/} | sed -E 's/R1/R2/' | sed -E 's/val_1/val_2/')
    echo $barcode","$dataset","$file1","$file2 >> ${META_OUT}
done
# END MAKE META FILE

# BEGIN MAKE CONFIGURATION FILE
echo "
# Needs to be run from wgbs-load/code/ directory!
# Needs to know about reference genome and index
reference = ../../../../reference/GENCODE/h38_no_alt.fa
index_dir = ../../../../reference/GENCODE/gembs-index/

sequence_dir = 01-fastq-trimmed
bam_dir = 02-mapping
bcf_dir = 03-calls
extract_dir = 04-extract
report_dir = 05-report

# Large memory footprint, less so on CPUs
memory = 80G
cores = 8
keep_logs = True

[calling]
right_trim = 0,0
left_trim = 0,0

[extract]
make_snps = False
make_cpg = False
make_non_cpg = False
make_bedmethyl = False" > ${CONF_OUT}
# END MAKE CONFIGURATION FILE


# GEMBS PREPARATION AND CONSOLE OUTPUT
gemBS prepare -c ${CONF_OUT} -t ${META_OUT}
echo "gemBS will run the following commands:"
gemBS --dry-run run
# END GEMBS PREP/OUTPUT


# MAPPING
parallel -S ${RUN_SERVERS} --joblog ${DATE_STR}-map.log --nonall --workdir . gemBS map
# END MAPPING


# Calling
parallel -S ${RUN_SERVERS} --joblog ${DATE_STR}-call.log --nonall --workdir . gemBS call
gemBS report
# END CALLING

# EXTRACT

mkdir -p -v "${EXTRACT_PATH}"
if [ -z "$(ls -A ${EXTRACT_PATH})" ] ;
then
    echo "${EXTRACT_PATH} is empty, extracting methylation now."

    # bcfs we want are in format 123.bcf
    ls -1 ${CALLS_PATH}???.bcf > INPUT
    # Pipes needed since there are '/' in the two variables
    ls -1 ${CALLS_PATH}???.bcf | sed -E 's/bcf/tsv/' | sed "s|$CALLS_PATH|$EXTRACT_PATH|g" > OUTPUT
  
    parallel --link --workdir . --joblog ${DATE_STR}-extract.log \
        ./extract.sh {1} {2} :::: INPUT :::: OUTPUT
else
   echo "${EXTRACT_PATH} is full (no need to extract methylation); delete if you need to!"
fi

cd ${MY_HOME}
done 