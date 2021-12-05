input_dir="../data/"
files="$(find ${input_dir} -type f -name "*[0-9][0-9][0-9].bcf")"

for infile in ${files}; do
    outfile=$(echo ${infile} | sed -E 's/03-calls/04-extract/' | sed -E 's/bcf/tsv/')

    # Get the file size (should be around 250 MB if it ran correctly)
    if [ -f "${outfile}" ]; then
        size=$(du -h ${outfile} | cut -f1 | sed 's/M//')
    else
        size="0"
    fi

    echo ${size}
done