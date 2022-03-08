# 02-call-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
    library(fdrtool)
    library(biomaRt)
})

# Argument parsing -----------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--idir", default= "./chr18/smooth-150-PCs-2/", help="directory with models.RData")
parser$add_argument("--alpha", default= 0.001, help="local false discovery rate cutoff")
args <- parser$parse_args()

load(file.path(args$idir, "models.RData"))
CHR <- test.cohort$chr[1]

# Get genes from ENSEMBL ------------------------------------------------------
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('chromosome_name','start_position','end_position', 
                    'hgnc_symbol', 'transcript_start', 'transcript_end'), 
                    filters = 'chromosome_name', values = str_remove(CHR, "chr"), mart = ensembl)
genes <- genes %>% dplyr::rename(chr = 'chromosome_name', start = 'start_position', end = 'end_position', gene.name = 'hgnc_symbol')
genes$chr = paste0("chr", genes$chr)

# HELPER FUNCTIONS -----------------------------------------------------------
.test_to_df <- function(test.cohort){
    # test.cohort (output from DSS) is a two-class object
    # grab the columns we need and output a cleaned up data frame
    df <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos, 
                stat = test.cohort$stat) %>% drop_na() 
    return(df)
}

.correct_pvals <- function(df){
    # given the processed test.cohort
    zz <- (df$stat - mean(df$stat)) / sd(df$stat)
    out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")

    # This one is how DSS calls DMRs...
    # Update: DSS uses areaStat so it's harder to hijack than I thought
    df$pvals <- out$pval

    # The rest are additional columns
    df$pval <- out$pval
    df$lfdr <- out$lfdr 
    df$qval <- out$qval
    return(df)
}

.callDMRs <- function(df, alpha=args$alpha){
    # Wrapper around DSS:callDMR that only allows alpha control
    class(df)[2] <- "DMLtest.multiFactor"
    dmrs <- callDMR(df, p.threshold=alpha, dis.merge = 2500, pct.sig = 0.5, minCG = 5) 
    return(dmrs)
}

.format_locus <- function(seq, a, b){
    # Given seq, a, b output chr6:12-15
    out <- paste0(seq, ":", as.character(a), "-", as.character(b))
    return(out)
}


.tally_sig_CpGs <- function(df, start, end, alpha=args$alpha){
    # Helper
    df %>% 
        filter(pos >= start) %>% 
        filter(pos <= end) %>%
        filter(pvals <= alpha) %>% 
        nrow() %>% return()
}

# Run functions to preapre data
df <- .test_to_df(test.cohort) %>% .correct_pvals()
dmrs <- df %>% .callDMRs()

# TODO: Vectorize
dmrs$nSigCG <- 0
for (i in 1:nrow(dmrs)){
    dmrs$nSigCG[i] <- .tally_sig_CpGs(df, dmrs$start[i], dmrs$end[i])
}


#--> How many are postivie/negative
.tally_signs <- function(df = diffs.df, start, end){
    signs <- df %>% filter(pos >= start) %>% filter(pos <= end) %>% pull(pos) %>% sign()
    nHyperMeth <- sum(signs == 1)
    return(nHyperMeth)
}


# TODO: Functionalize this---messy right now
window <- 10000
# Nearest gene
dmrs.gr <- GRanges(dmrs) 
genes.gr <- GRanges(genes)
overlap <- findOverlaps(dmrs.gr + window, genes.gr)

# Subset based on if they overlap
genes.gr.2 <- genes.gr[subjectHits(overlap), ]
dmrs.gr.2 <- dmrs.gr[queryHits(overlap), ]

# Compute distance between DMR and 10kb gene
dd <- as.data.frame(distanceToNearest(dmrs.gr.2, genes.gr.2))$distance

# DMR locus, distance, gene.locus, gene.name
dmrs.df <- as.data.frame(dmrs.gr.2, row.names = NULL) %>%
            transmute(dmr.locus = .format_locus(seqnames, start, end), length, nCG, areaStat, nSigCG)
dmrs.df$distance <- dd
ann.df <- as.data.frame(genes.gr.2, row.names = NULL) %>%
            transmute(gene.locus = .format_locus(seqnames, start, end), gene.name)
comb.df <- cbind(dmrs.df, ann.df)

write_csv(comb.df, file = file.path(args$idir, "closest-genes.csv"))

# Join the tally with the others
tmp <- comb.df %>%
    dplyr::select(dmr.locus, length, nCG, nSigCG, areaStat, distance, gene.name) %>%
    unique() %>% group_by(dmr.locus) %>% 
    dplyr::filter(gene.name != "") %>%
    summarize(dmr.locus, length, nCG, nSigCG, geneNames = paste0(gene.name, collapse = "; ")) %>% 
    dplyr::rename(dmrLocus = "dmr.locus") %>%
    unique() %>% arrange(-nCG)

write_csv(tmp, file = file.path(args$idir, "called-dmrs.csv"))


# Some plots
for (v in c("pval", "lfdr", "qval")){
    p <- df %>%
        ggplot(aes_string(x = v, y = "..density..")) +
        geom_histogram(bins = 50) +
        wiscR::light_theme()
    file = file.path(args$idir, paste0(v, ".png"))
    wiscR::save_plot(p, file)
}