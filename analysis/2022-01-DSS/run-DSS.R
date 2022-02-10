# run-DSS.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
library(data.table)
library(tidyverse)
library(argparse)
library(DSS)


parser <- ArgumentParser()
parser$add_argument("--chr", default= "chr6", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "NK,CD4T")
parser$add_argument("--num_pcs", default= "2", help = 'number of principal components (precomputed)')
args <- parser$parse_args()

# Setup output directory, will save called regions, figures?
odir <- paste0("result-", args$chr)
dir.create(odir, showWarning = FALSE)


pc.df <- read_csv(file.path("../../data/prin-comps-array-samples/", paste0(args$chr , ".csv")), 
                  col_types = cols())


ss.df <- read_csv("../../data/meta/phenos-cleaned.csv", col_types = cols())
DT <- fread(file.path("../../data/cov-meth/", paste0(args$chr, ".tsv")))


make_tall_matrix <- function(DT, value.var) {
    # Turn to wide format with samples indexing columns
    # and position / locus indexing row
    out <- dcast(DT, pos~sample, value.var=value.var)
    #rownames(out) <- out$pos
    #out[ ,c("pos") := NULL]
    
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


valid_samples <- intersect(intersect(colnames(M), ss.df$sample), pc.df$sample)

filt.df <- ss.df %>%
            inner_join(pc.df, by = "sample") %>% 
            dplyr::filter(sample %in% valid_samples) 

# pack years smoking needs to be not na
filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

drop_cols <- as.character(setdiff(colnames(M), valid_samples))[-1] # dont remove pos
M[, (drop_cols):=NULL]
Cov[, (drop_cols):=NULL]

print(length(valid_samples))

load.samples <- filt.df %>% filter(cohort == "AD") %>% pull(sample) %>% as.character 
ctrl.samples <- filt.df %>% filter(cohort == "CONTROL") %>% pull(sample) %>% as.character

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


print("Imputing values")
impute_by_group(M, load.samples)
impute_by_group(M, ctrl.samples)
impute_by_group(Cov, load.samples)
impute_by_group(Cov, ctrl.samples)

head(M)

M <- drop_null_positions(M, 0)
Cov <- drop_null_positions(Cov, 0)

# Order needs to be correct!!!
valid_pos <- intersect(M$pos, Cov$pos)
M <- M[pos %in% valid_pos , c("pos", filt.df$sample), with = FALSE]
Cov <- Cov[pos %in% valid_pos, c("pos", filt.df$sample), with = FALSE]

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(DT$chr[1], nrow(M)), pos = M$pos,
            M = as.matrix(M[ , -c("pos"), with=FALSE]), 
            Cov = as.matrix(Cov[, -c("pos"), with=FALSE]), 
            sampleNames = names(M)[-1])


# In[16]:


# Print this check in the future!!!!
all( filt.df$sample == colnames(bs) )


#TODO: formula from input flags
dml.fit <- DMLfit.multiFactor(bs, design = filt.df, smoothing = TRUE, smoothing.span = 200, 
            formula = ~cohort + PC1 + PC2 + Gran + CD8T + CD4T + NK + Bcell + bmi + age + sex)


# In[ ]:


colnames(dml.fit$X)
test.cohort <- DMLtest.multiFactor(dml.fit, coef = 2)


# In[ ]:


methylation <- M / Cov
methylation$pos <- M$pos
fwrite(methylation, file.path(odir, "methylation.csv"))


# In[ ]:


save(list = c("test.cohort"), file= file.path(odir, "test-values.RData"))
fwrite(M, file.path(odir, "M.csv"))
fwrite(Cov, file.path(odir, "Cov.csv"))
