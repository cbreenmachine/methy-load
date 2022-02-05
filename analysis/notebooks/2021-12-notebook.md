# December 2021 Notebook

One of my goals is to document everything better, espeically the history of this project. The format may change, but important things--keeping track of medium-major code changes, figuring out code idioms, downloading data a certain way--are all things worth documenting.

## December 15, 2021


## December 22, 2021

1. `trim_galore` only operating on directory
    - Problem: I had it set up so that `trim_galore` would run if a the `02-fastq-trimmed` directory did not exist or was empty. So, if four of five files were trimmed successfully, and I wanted to finish the fifith, I would delete the whole directory and re-run on all five. This is obviously inefficienct.
    - Solution: `runall.sh` now checks for *files* (based on `trim_galore` naming conventions) rather than *directories*.
2. `gemBS` freezing
    - Problem: when calling `gemBS map` in the `runall.sh` driver script, it would routinely freeze for 12 to 24 hours. It may be because of memory allocation, but I think the more likely problemn is that there are jobs that got "stuck" when `gemBS` was incompletely finished. 
    - Solution: added a `rm -rf .gemBS` line to `runall.sh` *before* the configuration file gets made. Seems to be working so far.
3. `munge.R`
    - Problem: hopefully none
