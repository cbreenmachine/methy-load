library(data.table)
library(dmrseq)
library(tidyverse)
library(annotatr)
library(dmrseq)
library(BiocParallel)
library(argparse)
register(MulticoreParam(6))

parser <- ArgumentParser()


parser$add_argument("--male-only", action="store_true", help="Filter out females")
parser$add_argument("--female-only", action="store_true", help="Filter out males")
parser$add_argument("-odir", default=paste0("figs-", Sys.Date(), "/"))
parser$add_argument("-chrom", default="22", help="Which chromosome to run analysis on.")
parser$add_argument("-covariates", default="CD8T,CD4T,NK,Bcell,Mono,Gran,age,bmi,sex,PC1,PC2")
parser$add_argument("-pheno.file", default="../../data/meta/phenos-cleaned.csv")
parser$add_argument("-mc.dir", default="../../data/cov_meth/")
parser$add_argument("-pc.dir", default="../../data/prin_comps/")
args <- parser$parse_args()
dir.create(args$odir, showWarning = FALSE)

# Derived variables
pc.file <- paste0(args$pc.dir, "PC_chr", args$chrom, ".tsv")
input.file <- paste0(args$mc.dir, "chr", args$chrom, "_cov_meth.tsv")
covariates <- as.vector(str_split(args$covariates, ",")[[1]])

# Put this into functions... / packages
# b3galt4 chr 6
# zadh2 chr 18

df <- read_csv(args$pheno.file) 
df <- df %>% inner_join(read_tsv(pc.file), by = "sample") %>% mutate(sample = as.character(sample))


block.covariate <- "sex"
# Handle sex filtering
if (args$male_only || args$chrom == "Y"){
    df <- df %>% dplyr::filter(sex == "M")
    covariates <- setdiff(covariates, "sex")
    block.covariate <- NULL
}

if (args$female_only){
    df <- df %>% dplyr::filter(sex == "F")
    covariates <- setdiff(covariates, "sex")
    block.covariate <- NULL
}


DT <- fread(input.file)
M <- dcast(DT, pos~sample, value.var="methylated")
C <- dcast(DT, pos~sample, value.var="coverage")

valid_samples <- intersect(unique(DT$sample), df$sample)

filt.df <- df %>%   dplyr::filter(sample %in% valid_samples) %>% 
                    column_to_rownames(var = "sample")

# pack years smoking needs to be not na
filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

drop_cols <- as.character(setdiff(colnames(M), valid_samples))[-1] # dont remove pos
M[, (drop_cols):=NULL]
C[, (drop_cols):=NULL]

# TODO: imputation...
M[is.na(M)] <- 0
C[is.na(C)] <- 0

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M)), pos = M$pos,
            M = as.matrix(M[ , -"pos", with=FALSE]), 
            Cov = as.matrix(C[, -"pos", with=FALSE]), 
            sampleNames = names(M)[-1])

pData(bs) <- filt.df %>%
                dplyr::select(c(cohort, machine, pool,
                                CD8T, CD4T, NK, Bcell, Mono, Gran, 
                                bmi, age, pack_years))

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)
sample.idx <- which(pData(bs)$cohort %in% c("AD", "CONTROL"))

bs.filtered <- bs[loci.idx, sample.idx]

regions <- dmrseq(bs = bs.filtered, cutoff = 0.01, testCovariate = "cohort", 
                    adjustCovariate = covariates, 
                    bpSpan = 50, minInSpan = 10)
# Aught to block on sex, 

save(regions, file = paste0(odir, "Y-regions.RData"))

annoTrack <- getAnnot("hg38")


plot_region <- function(regions, ix){
    ofile <- paste0(args$odir, "regions-", ix, ".png")
    png(ofile)

    plotDMRs(bs.filtered, regions=regions[ix,], 
         testCovariate="cohort",
         annoTrack=annoTrack)
    dev.off()
}

for (i in 1:100){
    plot_region(regions, i)
}