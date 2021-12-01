# munge.R
# converts between commonly used datatypes
# First iteration (November 30, 2021) goes from TSV to one big matrix

library(data.table)
library(tidyverse)
library(wiscR)
library(argparse)


ROOT_DIR <- "../data/"
samples.df <- read.table(file.path(ROOT_DIR, "meta/meta-data.tsv"), header=TRUE)
my_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].tsv", recursive = TRUE, full=TRUE)
print(file.info(my_files)$size)

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

# Need do this by methylated and unmethylated separate
DT <- rbindlist(lapply(my_files, read_tsv))

chr_range <-  paste0("chr", c(1:22, "X", "Y"))
for (z in chr_range){
    print(z)
    DT.tmp <- DT[which(chr == z), ]
    #DT.tmp[ ,chr:=NULL]
    fwrite(DT.tmp, paste0(z,"_meth_cov.tsv"), sep = "\t")


    # No need to subset, dcast will drop columns that aren't mentioned
    #dcast(DT.tmp, sample~pos, value.var="methylated") %>% fwrite(paste0(z,"_methylated.tsv"), sep = "\t")
    #dcast(DT.tmp, sample~pos, value.var="unmethylated") %>% fwrite(paste0(z,"_unmethylated.tsv"), sep = "\t")

}


