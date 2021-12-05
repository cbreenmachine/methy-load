#!/bin/bash
# extract_crawler.sh
# Given a parent directory (e.g. ../data/) it finds all 
input_dir="../data/"
files="$(find ${input_dir} -type f -name "*[0-9][0-9][0-9].bcf")"

> CRAWLER_INPUT
> CRAWLER_OUTPUT

for infile in ${files}; do
    outfile=$(echo ${infile} | sed -E 's/03-calls/04-extract/' | sed -E 's/bcf/tsv/')

    # Get the file size (should be around 250 MB if it ran correctly)
    if [ -f "${outfile}" ]; then
        size=$(du -h ${outfile} | cut -f1 | sed 's/M//')
    else
        size="0"
    fi
    # Test file size, add to queue if its small (or non-existent)
    if [ "$size" -gt "250"  ]; then
        echo "${outfile} already exists."
    else
        echo "${infile}" >> CRAWLER_INPUT
        echo "${outfile}" >> CRAWLER_OUTPUT
        new_dir=$(dirname "${outfile}")
        mkdir -vp "${new_dir}"
    fi
done

echo ""
echo "Will extract methylation on the following files:"
cat CRAWLER_INPUT

parallel --link --workdir . --joblog crawler.log \
    ./extract.sh {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT
