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
RUN_SERVERS="nebula-2,nebula-3,nebula-4"
DATE_STR=$(date +"%y-%m-%d")

# SETUP PATHS
#TODO: Make this accept from getopts
ROOT_DIR="../data/2021-11-01-illinois-batch1/"
POOL_SUB_DIRS="$(echo group{1..1}/)" # be sure to have trailing slash

for pool in ${POOL_SUB_DIRS}
do
    echo "Found these input files in ${ROOT_DIR}${pool}00-fastq/"
    ls ${ROOT_DIR}${pool}"00-fastq/"

    echo "Creating directory"
    mkdir -p -v ${ROOT_DIR}${pool}"01-fastq-trimmed/"
done

#TODO --joblog to parallel commands
echo "To continue running, hit y"
read -n 1 k <&1
if [[ $k = n ]] ; then
  echo "Re-configure directories or delete files as needed!"
elif [[ $k = y ]] ; then
  echo "continuing"
fi


for pool in $POOL_SUB_DIRS;
do
echo "Working on ${pool}"
meta_out="./$(echo ${pool%/}).meta.csv" # name of meta file used in gemBS
conf_out="./$(echo ${pool%/}).conf" # name of congiguration file used by gemBS
fastq_path=${ROOT_DIR}${pool}"00-fastq/"
fastq_trimmed_path=${ROOT_DIR}${pool}"01-fastq-trimmed/"

mkdir -p ${fastq_trimmed_path}

# TODO: Option to delete all from command line
# ./runall.sh --from_scratch 
# deletes all intermediate files

# trim files (conditional on them not being trimmed yet)
# Consider piping file names to TMP and then using parallel afterwords
ls -1 ${fastq_path}*R1*.fq | uniq > LEFT
ls -1 ${fastq_path}*R2*.fq | uniq > RIGHT

# 1. After trimming, rename files, then move them to appropriate dir
# 2. Make conditional ()
if [ -z "$(ls -A ${fastq_trimmed_path})" ]
then
   echo "${fastq_trimmed_path} is empty, trimming files now."
    # --link creates a mapping between the lines in LEFT and lines in RIGHT 
    # (one-to-one instead of pairwise combinations)
    # the fourth ':' means cat LEFT and RIGHT (don't treat as variable/expansion)
    parallel --link -S ${RUN_SERVERS} --workdir . \
        trim_galore --phred33 --cores 6 --output_dir ${fastq_trimmed_path} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT

    # (DEPRECATED) Move outputs of trim-galore into the correct folder *.fq *.html *.txt
    #for f1 in *.txt *.fq # add other extensions if doing fastqc
    #do
        # Remove adapter sequence in middle of file, also remove lane info from file name
        # What the hell is happening with these paths....?
     #   f2=$(echo ${f1} | sed -e 's/[-ACTG]//g' | sed -e 's/__L00M//g')
      #  mv "${f1}" "${fastq_trimmed_path}${f2}"
    #done

else
   echo "${fastq_trimmed_path} is full (no need to trim/fastqc files); delete if you need to re-trim!"
fi
# END TRIMMING SECTION


# BEGIN MAKE META FILE
echo "barcode,dataset,end1,end2" > ${meta_out}
for f in "${fastq_trimmed_path}"*val_1.fq
do
    # f="../../data/something.csv"
    # --> ${f##*/} is something.csv
    barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
    dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq//' | sed -E 's/_R1_//')
    file1=$(echo ${f##*/})
    file2=$(echo ${f##*/} | sed -E 's/R1/R2/' | sed -E 's/val_1/val_2/')
    echo $barcode","$dataset","$file1","$file2 >> ${meta_out}
done


# make configuration file
echo "
# Needs to be run from wgbs-load/code/ directory!
# Needs to know about reference genome and index
reference = ../../reference/GENCODE/h38_no_alt.fa
index_dir = ../../reference/GENCODE/gembs-index/


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
parallel -S ${RUN_SERVERS} --joblog ${ROOT_DIR}${pool}call.log --nonall --workdir . gemBS call
# END CALLING

gemBS report
# Extraction
# Use awk skript here

done 