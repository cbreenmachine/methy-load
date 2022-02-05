library(data.table)
library(tidyverse)
library(fastglm)
library(parallel)
library(gmodels)
library(argparse)

options(readr.show_col_types = FALSE)

# write covariates to a file in odir
parser <- ArgumentParser(description='Methylation estimates')
parser$add_argument('--covariates', default ="mean_methylation,CD8T,CD4T,NK,Bcell,Mono,Gran", help='where to store output figures.')
parser$add_argument('--ofile', default ="./data/PC-cell-compositions.csv", help='where to store output data.')
parser$add_argument('--sample', action ="store_true", help='where to store output data.')
args <- parser$parse_args()

covariates <- c(str_split(args$covariates,",")[[1]])

DT <- fread("../../data/cov_meth/chr22_cov_meth.tsv")
tmp.1 <- read_csv("../../data/meta/phenos-cleaned.csv")
tmp.2 <- read_tsv("../../data/prin_comps/PC_chr22.tsv")
samples.df <- merge(tmp.1, tmp.2, by = "sample")

# Avoid the scenario where methylation == 0|1
pseudocount <- 1

DT$methylation <- (DT$methylated + pseudocount) / (DT$coverage + 1)
DT.wide <- dcast(DT, sample~pos, value.var="methylation")

# number of observations and predictors
n <- length(unique(DT$sample))
p <- length(unique(DT$pos))
num.nas <- colSums(is.na(DT.wide))
DT.wide[ , names(which(num.nas > 0)):=NULL]


filter_DT <- function(DT, ix=2){
    # grabs first column (sample) and then some locus' methylation
    out.df <- as.data.frame(DT[ ,.SD,.SDcols=c(1,ix)])
    names(out.df)[2] <- 'y'

    out.df <- merge(out.df, samples.df, by = "sample")
    return(out.df)
}

fit_model <- function(input.df){
    # methylation point estimate only...
    x <- model.matrix( reformulate(covariates), data = input.df)
    y <- input.df$y

    glm.fitted <- fastglm::fastglm(y = y, x = x, family = quasibinomial(link = "logit"))
    y.out <- glm.fitted$fitted
    names(y.out) <- input.df$sample
    return(y.out)
}


# Computationally intensive part (first col is pos)
# If estimating the methylation percentage, then quasibinomial
# 3.6 minutes (my time) for 10,000

p <- 2:ncol(DT.wide)

if (args$sample){
    print("Sampling 10000 loci")
    set.seed(919)
    p <- sample(p, 10000, replace = FALSE)
}

system.time(
    out <- mclapply( p, function(x) (fit_model(filter_DT(DT.wide, x))), mc.cores = 8)
)

X <- do.call(cbind, out)
rownames(X) <- names(out[[1]])

pca.big <- fast.prcomp(X, center = TRUE, scale = FALSE) 
df <- as.data.frame(pca.big$x) %>% rownames_to_column(var = "sample")
write_csv(df, args$ofile)


