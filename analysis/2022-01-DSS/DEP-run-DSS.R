# run-DSS.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
# 

#TODO: Add to munging script that we want a wide M matrix and wide Cov matrix instead of one big thing...

library(data.table)
library(tidyverse)
library(argparse)
library(wiscR)
library(DSS)
library(annotatr)


parser <- ArgumentParser()
parser$add_argument("--chr", default= "chr6", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "")
parser$add_argument("--num_pcs", default= 2, help = 'number of principal components (precomputed)')
args <- parser$parse_args()

# Setup output directory, will save called regions, figures?
odir <- paste0("result-", args$chr)
dir.create(odir, showWarning = FALSE)


pc.df <- read_csv( file.path("../../data/prin_comps/", paste0(args$chr , ".tsv")) )


DT <- fread(args$ifile)
ss.df <- read_csv("../../data/meta/phenos-cleaned.csv")

M <- dcast(DT, pos~sample, value.var="methylated")
Cov <- dcast(DT, pos~sample, value.var="coverage")

valid_samples <- intersect(unique(DT$sample), ss.df$sample)

filt.df <- ss.df %>%
            inner_join(pc.df, by = "sample") %>% 
            dplyr::filter(sample %in% valid_samples) 

# pack years smoking needs to be not na
filt.df$pack_years[is.na(filt.df$pack_years)] <- 0

drop_cols <- as.character(setdiff(colnames(M), valid_samples))[-1] # dont remove pos
M[, (drop_cols):=NULL]
Cov[, (drop_cols):=NULL]

# TODO: imputation...
threshold <- floor( (ncol(M) - 1) / 2 )
threshold <- 0
drop.ix <- ( rowSums(is.na(M)) > threshold)
M <- M[!drop.ix, ]
Cov <- Cov[!drop.ix, ]
table(drop.ix)



# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(DT$chr[1], nrow(M)), pos = M$pos,
            M = as.matrix(M[ , -c("pos"), with=FALSE]), 
            Cov = as.matrix(Cov[, -c("pos"), with=FALSE]), 
            sampleNames = names(M)[-1])

dml.fit <- DMLfit.multiFactor(bs, design = filt.df, smoothing = TRUE, smoothing.span = 100, 
        formula = ~cohort + machine + CD8T + CD4T + NK + Bcell + Mono + bmi + age + sex)
colnames(dml.fit$X)
test.cohort <- DMLtest.multiFactor(dml.fit, coef = 2)

#save(list = c("dml.fit", "test.cohort"), file = file.path(args$odir, "fits.RData"))

# look at top sites
dmrs <- callDMR(test.cohort, p.threshold=0.01)



# NOT WORKING from here down

annoTrack <- getAnnot("hg38")




plot_region <- function(ix){

    start <- dmrs$start[ix]
    stop <- dmrs$end[ix]

    DT.sub <- DT[start < pos]
    DT.sub <- DT.sub[pos < stop]

    DT.sub$methylation <- DT.sub$methylated / DT.sub$coverage
    DT.wide <- dcast(DT.sub, sample~pos, value.var="methylation")

    p <- DT.sub %>%
        inner_join(filt.df, by = "sample") %>%
        ggplot(aes(x = pos, y = methylation, group = sample, color = cohort)) +
        geom_smooth(se = F) +
        geom_point(alpha = 0.5, size = 3) +
        wiscR::light_theme() +
        labs(caption = paste0("Length: ", dmrs$length[ix], "\t",
                              "nCG: ", dmrs$nCG[ix], "\t",
                              "areaStat:", round(dmrs$areaStat[ix], 2)))
    wiscR::save_plot(p, file.path(args$odir, paste0(start, ".png")))

    p <- DT.sub %>%
        inner_join(filt.df, by = "sample") %>%
        ggplot(aes(x = pos, y = methylation, group = cohort, color = cohort)) +
        geom_smooth(se = F) +
        geom_point(alpha = 0.5, size = 3) +
        wiscR::light_theme() +
        labs(caption = paste0("Length: ", dmrs$length[ix], "\t",
                              "nCG: ", dmrs$nCG[ix], "\t",
                              "areaStat:", round(dmrs$areaStat[ix], 2)))
    wiscR::save_plot(p, file.path(args$odir, paste0(start, "-group.png")))
}

for (ii in 1:nrow(dmrs)){
    plot_region(ii)
}

save(list = c("dmrs"), file = file.path(args$odir, "dmrs.RData"))

# TODO: 
a <- subsetByOverlaps()