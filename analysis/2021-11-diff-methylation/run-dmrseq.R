library(data.table)
library(dmrseq)
library(tidyverse)
library(annotatr)

# TODO: bring in other sample sheet for it's covariates.
# Adjust for PCs
# Put this into functions... / packages

# b3galt4 chr 6
# zadh2 chr 18
my_file <- "../../data/cov_meth/chr6_meth_cov.tsv"
DT <- fread(my_file)

tmp <- basename(my_file) %>% str_remove(".tsv") %>% str_split("_") 
chr <- tmp[[1]][1]

M <- dcast(DT, pos~sample, value.var="methylated")
C <- dcast(DT, pos~sample, value.var="coverage")

#M[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = 2:ncol(M)] / nrow(M)

# position vector should be sorted...it is
#tmp <- (M$pos[-1] - M$pos[1:(nrow(M)-1)]) > 0
#all(tmp)
samples.df <- read.csv("../../data/meta/archived/DHMRI_samplesheet_pilot.csv") 

valid_samples <- intersect(unique(DT$sample), samples.df$Alisch_ID)
samples_filt.df <- samples.df %>% dplyr::filter(Alisch_ID %in% valid_samples) %>%
                    dplyr::select(c(COHORT,AGE, GENDER, CD8T, CD4T, NK, Bcell, Mono, Gran))

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

pData(bs) <- samples_filt.df

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")==0) == 0)
sample.idx <- which(pData(bs)$COHORT %in% c("AD", "CONTROL"))

bs.filtered <- bs[loci.idx, sample.idx]

regions <- dmrseq(bs = bs.filtered, cutoff = 0.01, testCovariate = "COHORT")


annoTrack <- getAnnot("hg38")

png("tmp1.png")
plotDMRs(bs.filtered, regions=regions[1,], 
         testCovariate="COHORT",
         annoTrack=annoTrack)
dev.off()