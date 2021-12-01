
input_dir="../data/"
files="$(find ${input_dir} -type f -name "*[0-9][0-9][0-9].bcf")"

> CRAWLER_INPUT
> CRAWLER_OUTPUT

for infile in ${files}; do
    outfile=$(echo ${infile} | sed -E 's/03-calls/04-extract/' | sed -E 's/bcf/tsv/')
    num_lines=$(wc -l <"${outfile}")

    if [ -f "${outfile}" ] && [ "$num_lines" -gt "10000000"  ]; then
        echo "${outfile} already exists."
    else
        echo "${infile}" >> CRAWLER_INPUT
        echo "${outfile}" >> CRAWLER_OUTPUT

        new_dir=$(dirname "${outfile}")
        mkdir -p "${new_dir}"
    fi
done

parallel --link --workdir . --joblog crawler.log \
    ./extract.sh {1} {2} :::: CRAWLER_INPUT :::: CRAWLER_OUTPUT
