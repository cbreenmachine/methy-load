#!/bin/bash
# runall.sh
# conda environment: runs in wgbs (preferred for parallel cutadapt) or load
# Runs entire bioinformatics pipeline starting with paired-end data and ending with bcf files.
# Made parallel with GNU parallel

# Steps in pipeline
# 1. Trimming using _trim_galore and cutadapt
# 2. gemBS mapping
# 3. gemBS calling
# TODO 4:

# Practical considerations
# 1. gemBS is a high memory pipeline--stick to the higher memory servers
# 2. Benchhmarking to come.
# 3. Run on ../data/dummy/toy[1-2]/ first to make sure messages and logs are 
# behaving how they're supposed to. (structure is analogous to real data with 'pool' swapped for 'toy')

#TODO: Make default environment wgbs (after finishes running)
# CONSTANTS
RUN_SERVERS="nebula-4,nebula-7"
DATE_STR=$(date +"%y-%m-%d")

# SETUP PATHS
#TODO: Make this accept from getopts
ROOT_DIR="../data/2021-11-01-illinois-batch1/"
POOL_SUB_DIRS="$(echo toy{1..2}/)" # be sure to have trailing slash

for pool in ${POOL_SUB_DIRS}
do
    echo "Found these input files in ${ROOT_DIR}${pool}/00-fastq/"
    ls ${ROOT_DIR}${pool}"/00-fastq/"

    echo "Creating directory"
    mkdir -p -v ${ROOT_DIR}${pool}"/01-fastq-trimmed/"
done

#TODO --joblog to parallel commands
for sub_dir in $POOL_DIRS;
do
    ls ${ROOT_DIR}${sub_dir}

read -n 1 k <&1
if [[ $k = n ]] ;
then
exit ERRCODE "Re-configure directories or delete files as needed!"
elif [[ $k = y ]]; then
echo "continuing"

fi


for pool in $POOL_DIRS;
do
echo "Working on ${pool}"
meta_out="./${pool}-meta-tmp.csv" # name of meta file used in gemBS
conf_out="./${pool}-conf.conf" # name of congiguration file used by gemBS
fastq_path=${input_dir}${pool}"00-fastq/"
fastq_trimmed_path=${input_dir}"01-fastq-trimmed/"


#TODO: Make 
mkdir -p ${fastq_trimmed_path}

# Option to delete all from command line
# ./runall.sh --from_scratch 
# deletes all intermediate files

# trim files (conditional on them not being trimmed yet)
# Consider piping file names to TMP and then using parallel afterwords
ls -1 ${fastq_path}*R1*.gz | uniq > LEFT
ls -1 ${fastq_path}*R2*.gz | uniq > RIGHT

# 1. After trimming, rename files, then move them to appropriate dir
# 2. Make conditional ()
if [ -z "$(ls -A ${fastq_trimmed_path})" ]
then
   echo "${fastq_trimmed_path} is empty, trimming files now."
    # --link creates a mapping between the lines in LEFT and lines in RIGHT 
    # (one-to-one instead of pairwise combinations)
    # the fourth ':' means cat LEFT and RIGHT (don't treat as variable/expansion)
    parallel --joblog ./log/log.out --link -S ${RUN_SERVERS} \
        trim_galore --phred33 --cores 6 --output_dir ${fastq_trimmed_path} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT

    # (DEPRECATED) Move outputs of trim-galore into the correct folder *.fq *.html *.txt
    #for f1 in *.zip *.html *.txt *.fq
    #do
        # Remove adapter sequence in middle of file, also remove lane info from file name
    #    f2=$(echo ${f1} | sed -e 's/[-ACTG]//g' | sed -e 's/__L00M//g')
    #    mv "${f1}" "${fastq_trimmed_path}${f2}"
    #done

else
   echo "${fastq_trimmed_path} is full (no need to trim/fastqc files); delete if you need to re-trim!"
fi
# END TRIMMING SECTION


# BEGIN MAKE META FILE
echo "barcode,dataset,end1,end2" > ${pool}".meta.csv"
for f in "${fastq_trimmed_path}"*val_1.fq
do
    # f="../../data/something.csv"
    # --> ${f##*/} is something.csv
    barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
    dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq.gz//')
    file1=$(echo ${f##*/})
    file2=$(echo ${f##*/} | sed -E 's/R1_001/R2_001/' | sed -E 's/val_1/val_2/')
    echo $barcode","$dataset","$file1","$file2 >> ${meta_out}
done


# make configuration file
echo "
# Needs to be run from wgbs-load/code/ directory!
# Needs to know about reference genome and index
reference = ../../reference/GRCh38_latest_genomic.fna
index_dir = ../../reference/gembs-index/

sequence_dir = ${ROOT_DIR}${pool}01-fastq-trimmed
bam_dir = ${ROOT_DIR}${pool}02-mapping
bcf_dir = ${ROOT_DIR}${pool}03-calls
extract_dir = ${ROOT_DIR}${pool}04-extract
report_dir = ${ROOT_DIR}${pool}05-report

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
make_bedmethyl = False

" > ${conf_out}

# GEMBS PREPARATION AND CONSOLE OUTPUT
gemBS prepare -c ${conf_out} -t ${meta_out}
echo "gemBS will run the following commands:"
gemBS --dry-run run
# END GEMBS PREP/OUTPUT


# MAPPING
parallel -S ${RUN_SERVERS} --joblog ${ROOT_DIR}${pool}map.log --nonall --workdir . gemBS map
# END MAPPING


# Calling
parallel -S ${RUN_SERVERS} --joblog ${ROOT_DIR}${pool}call.log --nonall --workdir . gemBS map
# END CALLING


# Extraction
# Use awk skript here

done # END LOOP THROUGH SUB DIRECTORIES