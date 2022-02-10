library(data.table)
library(tidyverse)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--ifile", default= "../../cov-meth/chr22.tsv", help='Chromosome to run DSS on')
parser$add_argument("--filter", action = "store_true")
parser$add_argument("--filter.file", default= "../../data/meta/array-samples.csv")
parser$add_argument("--odir.M", default= "../../data/imputed-coverage/")
parser$add_argument("--odir.Cov", default= "../../data/imputed-methylated/")
args <- parser$parse_known_args()

dir.create(args$odir.M, showWarnings = FALSE)
dir.create(args$odir.Cov, showWarnings = FALSE)


make_tall_matrix <- function(DT, value.var) {
    # Turn to wide format with samples indexing columns
    # and position / locus indexing row
    out <- dcast(DT, pos~sample, value.var=value.var)
    return(out)
}


drop_null_positions <- function(X, p=0.5){
    # Drops positions in X with too many nulls,
    # X is from `make_tall_matrix` and p is a percentage from 0 to 1
    stopifnot(p >= 0, p <= 1)
    # Cut point as float
    cut <- (ncol(X) - 1) * p
    num_null <- rowSums(is.na(X))
    
    return(X[num_null <= cut, ])
}

get_na_mask <- function(X, col.names){
    mask <- is.na(X)
    cn <- colnames(mask)
    for (j in 1:ncol(mask)){
        if ((cn[j] %in% col.names)) {
            mask[, j] <- FALSE
        }
    }
    return(mask)
}

impute_by_group <- function(X, col.names, round.mean = TRUE){
    # The imputed value is the position mean
    position.means <- rowMeans(X[ , col.names, with = FALSE], na.rm = TRUE)
    
    # If 
    if (round.mean) {position.means <- round(position.means)}
    mask <- get_na_mask(X, col.names)
    
    # Imputation step
    # With data.tables, for loops much master than vectorization
    for (i in 1:nrow(X)){
        for (j in which(mask[i, ])){
            set(X, i, j, value = position.means[i])
        }
    }
}






filter.samples <- read_csv(args$filter.file) %>% pull(sample)
DT <- fread(args$ifile)

M <- make_tall_matrix(DT, "methylated") %>% drop_null_positions
Cov <- make_tall_matrix(DT, "coverage") %>% drop_null_positions

# Samples we have data on meet samples we want to process
valid.samples <- intersect(colnames(M), filter.samples)
print(length(valid.samples))

# Drop the ones we don't want
drop_cols <- as.character(setdiff(colnames(M), valid.samples))[-1] # dont remove pos
M[, (drop.cols):=NULL]
Cov[, (drop.cols):=NULL]

load.samples <- filt.df %>% filter(cohort == "AD") %>% pull(sample) %>% as.character 
ctrl.samples <- filt.df %>% filter(cohort == "CONTROL") %>% pull(sample) %>% as.character
















