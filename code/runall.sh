#!/bin/bash

# runall.sh


# SETUP PATHS
ROOT_DIR="../data/2021-11-01-illinois-batch1/"
POOL_DIRS="$(echo pool{1..4})"


#TODO: print times, auto print to log
#TODO: better path specification. while getopts with flags. Decide after batch file structure is settled


#TODO: Check paths. Creating 
echo "The following input fastq files are identified:"

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



for thing in $POOL_DIRS;
do
  echo $thing
done

meta_out="./meta-tmp.csv"

# Probably want to define variables here
fastq_path=${input_dir}"00-fastq-110-only/"
fastq_trimmed_path=${input_dir}"01-fastq-trimmed-110-only/"

mkdir -p ${fastq_trimmed_path}

# Option to delete all from command line
# ./runall.sh --from_scratch 
# deletes all intermediate files

# trim files (conditional on them not being trimmed yet)
# Consider piping file names to TMP and then using parallel afterwords
ls -1 ${fastq_path}*R1*.gz | uniq > LEFT
ls -1 ${fastq_path}*R2*.gz | uniq > RIGHT

# TODO: This is bottlenecked because cutadapt can only run on one core. Fine when processing multiple files, but annoying when just one.
# 1. After trimming, rename files, then move them to appropriate dir
# 2. Make conditional ()
if [ -z "$(ls -A ${fastq_trimmed_path})" ]
then
   echo "${fastq_trimmed_path} is empty, trimming files now."
    # NODEFILE tells GNU parallel which servers it can use
    # --link creates a mapping between the lines in LEFT and lines in RIGHT
    # which is one line == one line (as opposed to pairwise combos)
    # the fourth : means read LEFT and RIGHT (don't treat as variable/expansion)
    parallel --link \
        trim_galore --fastqc --phred33 \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT

    # Move outputs of trim-galore into the correct folder *.fq *.html *.txt
    for f1 in *.zip *.html *.txt *.fq
    do
        # Remove adapter sequence in middle of file, also remove lane info from file name
        f2=$(echo ${f1} | sed -e 's/[-ACTG]//g' | sed -e 's/__L00M//g')
        mv "${f1}" "${fastq_trimmed_path}${f2}"
    done

else
   echo "${fastq_trimmed_path} is full (no need to trim/fastqc files); delete if you need to re-trim."
fi
# END TRIMMING SECTION


# BEGIN MAKE META FILE
echo "barcode,dataset,end1,end2" > $meta_out
for f in "${fastq_trimmed_path}"*val_1.fq
do
    # f="../../data/something.csv"
    # --> ${f##*/} is something.csv
    barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
    dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq.gz//')
    file1=$(echo ${f##*/})
    file2=$(echo ${f##*/} | sed -E 's/R1_001/R2_001/' | sed -E 's/val_1/val_2/')
    echo $barcode","$dataset","$file1","$file2 >> $meta_out
done


# make configuration file
echo "
# Needs to know about reference genome and index
reference = ../../reference/GRCh38_latest_genomic.fna
index_dir = ../../reference/gembs-index/

sequence_dir = ${input_dir}01-fastq-trimmed
bam_dir = ${input_dir}02-mapping
bcf_dir = ${input_dir}03-calls
extract_dir = ${input_dir}04-extract
report_dir = ${input_dir}05-report

memory = 60G
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

" > conf.conf


# GEMBS PREPARATION AND CONSOLE OUTPUT
gemBS prepare -c conf.conf -t meta.csv
echo "gemBS will run the following commands:"
gemBS --dry-run run
# END GEMBS PREP/OUTPUT


# MAPPING
parallel -S nebula-2, nebula-5 --nonall --workdir . gemBS map
# END MAPPING


# Calling
#parallel -S nebula-3,nebula-5,nebula-7 --nonall ::: 'gemBS call'
gemBS call
# END CALLING


# Extraction
# Use awk skript here