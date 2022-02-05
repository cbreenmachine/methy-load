#! /bin/bash
# Coverage numbers for 
find ../../data/batch01/ -name "[0-9][0-9][0-9].bcf" > INPUT

function compute_stats() {
    # Input output files
    input=${1}
    output=$(basename ${1} .bcf)".stats"

    bcftools query --format '[%MC8]\n' --include 'CG="Y"' ${input} \
        | /s/bin/Rscript -e "summary(rowSums(read.csv(file(\"stdin\"))))" > ${output}
}

export -f compute_stats

parallel --jobs 8 --link compute_stats :::: INPUT
rm INPUT

echo "sample,min,q1,median,mean,q3,max" > coverage.csv

for ff in *.stats
do
    # Remove weird spaces from R output, replace with commas
    # strip first and last chracters from line
    line=$(cat ${ff} | tail -1 | sed 's/  */,/g' | tail -c +2 | head -c -2)

    # If the 
    if [ ! -z "${line}" ]
    then
        echo $(basename ${ff} .stats)","${line} >> coverage.csv
    fi
done
