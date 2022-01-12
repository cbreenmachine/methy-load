library(data.table)
library(dmrseq)
library(tidyverse)
library(annotatr)
library(BiocParallel)
library(argparse)
register(MulticoreParam(6))

#TODO: argparse, put multiple runs in bash script
#Allow filtering sex


odir <- paste0("figs-", Sys.Date(), "/")
dir.create(odir)

# Put this into functions... / packages
# b3galt4 chr 6
# zadh2 chr 18

df <- read_csv("../../data/meta/phenos-cleaned.csv") %>% 
        mutate(sample = as.character(sample)) %>%
        dplyr::filter(sex == "M")
my_file <- "../../data/cov_meth/chrY_cov_meth.tsv"
DT <- fread(my_file)

tmp <- basename(my_file) %>% str_remove(".tsv") %>% str_split("_") 
chr <- tmp[[1]][1]

M <- dcast(DT, pos~sample, value.var="methylated")
C <- dcast(DT, pos~sample, value.var="coverage")


valid_samples <- intersect(unique(DT$sample), df$sample)

filt.df <- df %>%   dplyr::filter(sample %in% valid_samples) %>% 
                    column_to_rownames(var = "sample")

filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

drop_cols <- as.character(setdiff(colnames(M), valid_samples))[-1] # dont remove pos

M[, (drop_cols):=NULL]
C[, (drop_cols):=NULL]

M[is.na(M)] <- 0
C[is.na(C)] <- 0

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(chr, nrow(M)), pos = M$pos,
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

regions <- dmrseq(bs = bs.filtered, cutoff = 0.01, testCovariate = "cohort")
save(regions, file = paste0(odir, "Y-regions.RData"))

annoTrack <- getAnnot("hg38")




plot_region <- function(regions, ix){
    ofile <- paste0(odir, "regions-", ix, ".png")
    png(ofile)

    plotDMRs(bs.filtered, regions=regions[ix,], 
         testCovariate="cohort",
         annoTrack=annoTrack)
    dev.off()
}

for (i in 1:100){
    plot_region(regions, i)
}