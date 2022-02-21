# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
    library(dmrseq)
})

parser <- ArgumentParser()
# Things that may will change
parser$add_argument("--chr", default= "chr6", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "Gran,CD8T,CD4T,NK,Bcell,bmi,age,sex")
parser$add_argument("--num_pcs", default= 2, help= 'Number of principal components to include in analysis')
parser$add_argument("--smoothing", default= 100, help= 'Width of smoothing window')
# Flags for methods
parser$add_argument("--DSS", action = "store_true", help = "Which method to use")
parser$add_argument("--DMRseq", action="store_true", help = "Which method to use")
# File paths
parser$add_argument("--samples_file", default= "../../data/meta/array-samples.csv", help="CSV file with samples to filter")
parser$add_argument("--pheno_file", default= "../../data/meta/phenos-cleaned.csv", help="CSV file with LOAD/control status")
parser$add_argument("--pc_dir", default= "../../data/prin-comps-array-samples/", help="CSV file with LOAD/control status")
args <- parser$parse_args()

covariates <- as.vector(unlist(str_split(args$covariates, ",")))
if (args$num_pcs > 0){
    covariates <- c(covariates, paste0("PC", 1:args$num_pcs))
}
dss.formula <- as.formula(paste(c("~cohort", covariates), collapse = "+"))


idir <- file.path("./", args$chr)
odir <- paste0(args$chr, "/", "smooth-", as.character(args$smoothing), "-PCs-", args$num_pcs, "/")
dir.create(odir, showWarnings = FALSE)

pc.df <- read_csv(file.path(args$pc_dir,  paste0(args$chr, ".csv")), col_types = cols())
ss.df <- read_csv(args$pheno_file, col_types = cols())

M <- fread(file.path(idir, "M.csv"))
Cov <- fread(file.path(idir, "Cov.csv"))

valid.samples <- intersect(intersect(colnames(M), ss.df$sample), pc.df$sample)

filt.df <- ss.df %>%
            inner_join(pc.df, by = "sample") %>% 
            dplyr::filter(sample %in% valid.samples) 

drop.cols <- as.character(setdiff(colnames(M), valid.samples))[-1] # dont remove pos

if (length(drop.cols) > 0){
    M[, (drop.cols):=NULL]
    Cov[, (drop.cols):=NULL]
}

#--> Handle missing values?
#filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

# Handle mismatching files
all(filt.df$sample == names(M)[-1])

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), pos = M$pos,
            M = as.matrix(M[ , -c("pos"), with=FALSE]), 
            Cov = as.matrix(Cov[, -c("pos"), with=FALSE]), 
            sampleNames = names(M)[-1])

all( filt.df$sample == colnames(bs) )


# Derive some parameters
smooth = TRUE
if (args$smoothing == 0){
    smooth = FALSE
}

if (args$DSS){

    dss.time <- system.time(
        dml.fit <- DMLfit.multiFactor(bs, design = filt.df, smoothing = smooth, 
        smoothing.span = args$smoothing, formula = dss.formula)
    )
    test.var <- colnames(dml.fit$X)[2]
    print(test.var)
    test.cohort <- DMLtest.multiFactor(dml.fit, coef = test.var)
}

if (args$DMRseq){
    #--> DMR seq version of analysis
    pData(bs) <- filt.df %>% dplyr::select(c("cohort", all_of(covariates)))

    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)
    sample.idx <- which(pData(bs)$cohort %in% c("AD", "CONTROL"))
    bs.filtered <- bs[loci.idx, sample.idx]

    register(MulticoreParam(4)) # Restrict because otherwise not enough memory...
    dmr.time <- system.time(
        regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, 
                        testCovariate = "cohort", 
                        minNumRegion = 5, 
                        adjustCovariate = covariates, 
                        bpSpan = args$smoothing, minInSpan = 10)
    )
}

#rm(M, Cov, bs, valid.samples, pc.df, drop.cols, loci.idx, sample.idx)

# Save the models
outname <- file.path(odir, "models.RData")
save(list = intersect(ls(), c("test.cohort", "regions", "dss.time", "dmr.time", "filt.df")), file = outname)