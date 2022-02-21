# run-DSS.R
#TODO: turn this into python--should be part of the python postprocess script...
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--idir", default = "../../data/cov-meth/", help = "where the ")
parser$add_argument("--chr", default= "chr6", help="Chromosome on which to impute")
parser$add_argument("--samples_file", default= "../../data/meta/array-samples.csv", help="CSV file with samples to filter")
parser$add_argument("--pheno_file", default= "../../data/meta/phenos-cleaned.csv", help="CSV file with LOAD/control status")
args <- parser$parse_args()

# Setup output directory, will save called regions, figures?
odir <- paste0(args$chr)
dir.create(odir, showWarning = FALSE)

# Read in data
keep.samples <- read_csv(args$samples_file, show_col_types = FALSE) %>% pull(sample) %>% as.character
ss.df <- read_csv(args$pheno_file, show_col_types = FALSE) %>% filter(sample %in% keep.samples)

#TODO: Change to csv
DT <- fread(file.path(args$idir, paste0(args$chr, ".tsv")))

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

M <- make_tall_matrix(DT, "methylated") %>% drop_null_positions
Cov <- make_tall_matrix(DT, "coverage") %>% drop_null_positions
print("Read in data and made tall")

valid.samples <- intersect(colnames(M), ss.df$sample)
filt.df <- ss.df %>% dplyr::filter(sample %in% valid.samples) 



drop_cols <- as.character(setdiff(colnames(M), valid.samples))[-1] # dont remove pos
M[, (drop_cols):=NULL]
Cov[, (drop_cols):=NULL]

print(length(valid.samples))

load.samples <- filt.df %>% filter(cohort == "AD") %>% pull(sample) %>% as.character 
ctrl.samples <- filt.df %>% filter(cohort == "CONTROL") %>% pull(sample) %>% as.character

get_na_mask <- function(X, col.names){
    # X is a data table with columns indexed by samples, rows by position
    # col.names specifies which columns we wish to impute
    # usually this is LOAD samples or control samples
    mask <- is.na(X) # TRUE --> will impute
    cn <- colnames(mask)
    for (j in 1:ncol(mask)){

        if ((cn[j] %in% col.names)) {
            mask[, j] <- FALSE # don't include others
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

print("Imputing values")
impute_by_group(M, load.samples)
impute_by_group(M, ctrl.samples)
impute_by_group(Cov, load.samples)
impute_by_group(Cov, ctrl.samples)

M <- drop_null_positions(M, 0)
Cov <- drop_null_positions(Cov, 0)

# Order needs to be correct!!!
valid.pos <- intersect(M$pos, Cov$pos)
M <- M[pos %in% valid.pos , c("pos", filt.df$sample), with = FALSE]
Cov <- Cov[pos %in% valid.pos, c("pos", filt.df$sample), with = FALSE]

# Print this check in the future!!!!
all(filt.df$sample == colnames(M)[-1])

methylation <- M / Cov
methylation$pos <- M$pos

# Write it all out
fwrite(methylation, file.path(odir, "methylation.csv"))
fwrite(M, file.path(odir, "M.csv"))
fwrite(Cov, file.path(odir, "Cov.csv"))
#END