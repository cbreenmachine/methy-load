# load-methy todo

Ways to improve code readability, function, and portability.

### Todo
- [ ] Bismark validation
  - [ ] Run on one sample in every group
- [ ] Link to bioinformatics style guide, etc.
- [ ] Remove delete permissions from raw files? See [here](https://ostechnix.com/prevent-files-folders-accidental-deletion-modification-linux/)

### In Progress
- [ ] Fix this in meta.csv files: 217,217_CGCAACTA-GAATCCGA_L00M001_,217_CGCAACTA-GAATCCGA_L00M_R1_001_val_1.fq,217_CGCAACTA-GAATCCGA_L00M_R2_001_val_2.fq
- [ ] Bismark reference genome on same NIH reference.
- [ ] Incorporate sandbox into `runall.sh`
- [ ] Console and logging output for `runall.sh`
- [ ] Tracking data delivery and notes in spreadsheet
- [ ] Protect raw files from accidental deletion

### Done ✓
- [x] Add "-" to log files in `runall.sh`
- [x] Incorporate `infer-methylation.awk` into pipeline. Or make separate?
- [x] cd to child directory and rum `gemBS` from there, or keep all configurations in one place?
- [x] Generate dated log files for trimming, mapping, calling
- [x] Automate pipeline and create working `runall.sh`
- [x] Parallelize trimming (across and within servers)