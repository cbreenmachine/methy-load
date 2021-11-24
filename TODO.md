# load-methy todo

Ways to improve code readability, function, and portability.

### Todo
- [ ] Experiment with differential analysis--DSS and 
  - [ ] Run on one sample in every group
- [ ] Link to bioinformatics style guide, etc.
- [ ] Remove delete permissions from raw files? See [here](https://ostechnix.com/prevent-files-folders-accidental-deletion-modification-linux/)

### In Progress
- [ ] Refactor R functions and put extract scripts into a module?
- [ ] Console and logging output for `runall.sh`
- [ ] Protect raw files from accidental deletion

### Done ✓
- [x] Incorporate sandbox into `runall.sh`
- [x] Bismark validation
- [x] PCA Plots
- [x] Add "-" to log files in `runall.sh`
- [x] Incorporate `infer-methylation.awk` into pipeline. Or make separate?
- [x] cd to child directory and rum `gemBS` from there, or keep all configurations in one place?
- [x] Generate dated log files for trimming, mapping, calling
- [x] Automate pipeline and create working `runall.sh`
- [x] Parallelize trimming (across and within servers)
- [x] Fix this in meta.csv files: 217,217_CGCAACTA-GAATCCGA_L00M001_,217_CGCAACTA-GAATCCGA_L00M_R1_001_val_1.fq,217_CGCAACTA-GAATCCGA_L00M_R2_001_val_2.fq
- [x] Bismark reference genome on same NIH reference.
- [x] Tracking data delivery and notes in spreadsheet