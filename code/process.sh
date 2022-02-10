

source configuration.txt


trim() {
  echo "Starting trimming process now"
  echo ${CONF_TXT}
}




prepare() {
  cd ${idir}
  # BEGIN MAKE META FILE
  echo "barcode,dataset,end1,end2" > ${META_OUT}
  for f in "${FASTQ_TRIMMED_PATH}"*val_1.fq
  do
      # f="../../data/something.csv"
      # --> ${f##*/} is something.csv
      barcode=$(echo ${f##*/} | cut -d "_" -f 1 | sed -E 's/R/s/')
    # dataset=$(echo ${f##*/} | sed -E 's/val_[1-2].fq//' | sed -E 's/_R1_//')
      file1=$(echo ${f##*/})
      file2=$(echo ${f##*/} | sed -E 's/R1/R2/' | sed -E 's/val_1/val_2/')
      echo "${barcode},${barcode},${file1},${file2}" >> ${META_OUT}
  done
  # END MAKE META FILE
  echo "Preparing gemBS, should only take a second"
  echo "${CONF_TXT}" > ${CONF_OUT}
  # GEMBS PREPARATION AND CONSOLE OUTPUT
  rm -rf .gemBS # clears some hanging errors
  gemBS prepare -c ${CONF_OUT} -t ${META_OUT}
  echo "gemBS will run the following commands:"
  gemBS --dry-run run
  cd ${home}
}

map() {
  # MAPPING
  cd ${idir}
  parallel -S mastodon-3 --joblog ${date}-map.log --nonall --workdir . gemBS map
  cd ${home}
  # END MAPPING
}

call(){
  # Calling
  parallel -S ${RUN_SERVERS} --joblog ${date}-call.log --nonall --workdir . gemBS call
  gemBS report
  # END CALLING
}

# if invoked as a script rather than sourced, call function named on argv via the below;
# note that this must be the first operation other than a function definition
# for $_ to successfully distinguish between sourcing and invocation:
#[[ $_ != $0 ]] && return

# make sure we actually *did* get passed a valid function name
if declare -f "$1" >/dev/null 2>&1; then
  # invoke that function, passing arguments through
  "$@" # same as "$1" "$2" "$3" ... for full argument list
else
  echo "Function $1 not recognized" >&2
  exit 1
fi
