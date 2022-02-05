# munge.R
# converts between commonly used datatypes
# First iteration (November 30, 2021) goes from TSV to one big matrix

library(data.table)
library(tidyverse)
library(wiscR)
library(argparse)
library(parallel)

ROOT_DIR <- "../data/"
samples.df <- read.table(file.path(ROOT_DIR, "meta/meta-data.tsv"), header=TRUE)
my_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].tsv", recursive = TRUE, full=TRUE)
odir <- "../../data/2021-12-02-methylation-coverage/"

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--pca', action='store_true', help='load raw files and recompute PCA')
args <- parser$parse_args()

##### FUNCTIONS ####
read_tsv <- function(file_name){
    # Read using fread, estimate methylation
    sample_name <- str_remove(file_name, ".tsv") %>% basename()

    DT <- fread(file_name)
    DT[ ,coverage:=methylated + unmethylated]
    DT.filtered <- DT[which(coverage > 0), ] # do this without creating new object?

    DT.filtered$sample <- sample_name
    
    DT.filtered[, c("reference","context", "strand"):=NULL]  # remove columns
    DT.filtered$pos <- as.character(DT.filtered$pos)
    print(paste("Read in", sample_name))
    return(DT.filtered)
}
#[which(chr == chr_select), ]
#DT <- read_tsv(my_files[1:3])

# Need do this by methylated and unmethylated separate
chr_range <-  paste0("chr", c(1:22, "X", "Y"))

# This needs options passed to read_tsv...
DT <- rbindlist(mclapply(my_files[1:3], read_tsv))

mclapply(chr_range, function(z) fwrite(DT[which(chr == z), ], paste0(odir, z,"_meth_cov.tsv"), sep = "\t"), mc.cores = 8)