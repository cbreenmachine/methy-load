# load-methy todo

Ways to improve code readability, function, and portability.

### Todo

- [ ] Incorporate `infer-methylation.awk` into pipeline. Or make separate?
- [ ] Bismark validation
  - [ ] Run on one sample in every group
- [ ] Rename directories in `./data/2021-11-01-illinois-batch1` to `poolxx-group[1-4]`
- [ ] Link to bioinformatics style guide, etc.

### In Progress

- [ ] Console and logging output for `runall.sh`
- [ ] Tracking data delivery and notes in spreadsheet
- [ ] Protect raw files from accidental deletion

### Done ✓

- [x] Generate dated log files for trimming, mapping, calling
- [x] Automate pipeline and create working `runall.sh`
- [x] Parallelize trimming (across and within servers)