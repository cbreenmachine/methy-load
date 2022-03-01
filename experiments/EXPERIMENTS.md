
# experiments.md

This file is my attempt at keeping a lab notebook. I try to document what did and didn't work (code, science, etc.). It is ordered chronologically by month and then thematically. This may change over time.

# October 2021
## GNU Parallel, SSH
- In my `~/.bashrc` file, I've added a command `source activate load` to make the default conda environment the one used in the LOAD project.
- Added the private key (shared with all nebula servers) to my list of known hosts. Now I can ssh without password.

### GNU Parallel example commands

This one says run the same command on all servers listed after -S.
`parallel --nonall -S nebula-2,nebula-3 lscpu`

### Cutadapt multi-core
cutadapt shoul work with muti-core, but doesn't. Boils down to needing to reduce python version (3 needed) from header, but the header shebang chagned to a bash flavored shbang. So basically it can't figure out that you're using python 3, so doesn't set pigz (parallel gzip) to on.

# November 2021

For whole-genome data, I am also using gnu parallel. Here's an idiom that will be useful for more complicated pipes. 

1. Break your pipe into a function that works with one file.
`fastq_to_sorted() {
    prefix="$1"
    output_dir="$2"
    out="${output_dir}${tmp}"
    do_stuff ${prefix}.in.csv ${out}   
}`
2. Export the function (not sure exactly why this is needed--probably a subshell thing)
`export -f fastq_to_sorted`
3. Run parallel. As best as I can tell, `:::` is a pipe (but call to parallel happens in the beginning instead of the end). The fourth colon means "read this file line by line" (or whatever delimiter you specify) instead of interpreting it literally. `parallel` will make pairwise combinations, although since we're doing n numer of files by 1, don't need to worry about this.
`parallel fastq_to_sorted :::: PREFIX ::: "${output_dir}"`

## Multiple Servers with Conda Environments

Sub shell needs to know about conda. Try running the following.
`source ~/miniconda3/etc/profile.d/conda.sh && conda activate`

## GNU Parallel Errors

Running `parallel -S nebula-1,nebula-2 --nonall "gemBS map"` resulted in the following error:

`parallel: Error: -g has been retired. Use --group.
parallel: Error: -B has been retired. Use --bf.
parallel: Error: -T has been retired. Use --tty.
parallel: Error: -U has been retired. Use --er.
parallel: Error: -W has been retired. Use --wd.
parallel: Error: -Y has been retired. Use --shebang.`


## Library Prep Kit
Illinois uses [NEBNext Enzymatic Methyl-seq Kit](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit#Product%20Information). The technical note is found [here](https://www.neb.com/-/media/nebus/files/application-notes/technote_nebnext_enzymatic_methyl-seq.pdf?rev=015e017f782a4bc9b50b61f1b0f3c807&hash=E23FFED5E02C7A282E753B6B8F555A30). It is also downloaded and stored in `resources/`  as "TechNote_NEBNext...".

- [EM-seq libraries are directional](https://www.neb.com/faqs/2019/03/12/are-em-seq-libraries-directional-or-non-directional)
- 21 million CpGs common to all libraries in technical note.