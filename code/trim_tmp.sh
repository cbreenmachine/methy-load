
clean_trim_dir() {

  cd ${idir}
  cd ${FASTQ_TRIMMED_PATH}
  # remove anything that didn't finish
  rm *trimmed.fq
  unique_nums=$(find . -name "*.txt" | cut -c3-5) # get list like 235 236 237 ...

  for un in ${unique_nums}; 
  do 
    if [[ $(ls ${un}* | wc -l) -lt 4 ]]; then
        echo "Deleting files starting with ${un}"
        rm ${un}*
    fi
  done
  cd ${home}

}

compress_trim_dir() {
    # Check that extract worked?
    find "${FASTQ_TRIMMED_PATH}" \( -name '*.fq' \) -exec gzip --verbose --keep {} \;
}

trim() {
    cd ${idir}
    > LEFT # LEFT reads stored line-by-line in this file
    > RIGHT
    for ff in ${FASTQ_PATH}*R1*.fastq.gz; do 
        un=$(basename ${ff} | cut -c1-3) # grab the 123 from 00-fastq/123_stuff.fastq.gz

        # a successful trim_galore run produces 4 files starting with 123
        # if less than 4, needs trimming
        if [[ $(find "${FASTQ_TRIMMED_PATH}" -name "${un}*" | wc -l) -lt 4 ]] 
        then
            echo "${ff}" >> LEFT; 
            echo ${ff} | sed -e s/R1/R2/ | sed -e s/val_1/val_2/ >> RIGHT; 
        fi
    done

    # --link creates a mapping between the lines in LEFT and lines in RIGHT 
    # (one-to-one like python's zip instead of pairwise combinations)
    # the fourth ':' means cat LEFT and RIGHT (don't treat as variable/expansion)
    parallel --link -S ${RUN_SERVERS} --workdir . --joblog ${DATE_STR}-trim.log \
        trim_galore --phred33 --cores 6 --output_dir ${FASTQ_TRIMMED_PATH} \
        --dont_gzip --paired {1} {2} :::: LEFT :::: RIGHT
    rm LEFT RIGHT # cleanup

    cd ${home}
}