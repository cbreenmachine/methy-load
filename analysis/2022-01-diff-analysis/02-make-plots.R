# make-plots.R
# Plots tracks for three types of data
# 1. Newly found DMPs
# 2. Regions of interest
# 3. Concordance with other DMPs
suppressPackageStartupMessages({
    library(DSS)
    library(Gviz)
    library(tidyverse)
    library(data.table)
    library(ggsci)
    library(parallel)
    library(scales)
    library(wiscR)
    library(argparse)
    library(methanalyzeR)
    library(Homo.sapiens)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}) 

#PLOTS WISHLIST
# 1. Comprehensive plot of DMRs from DSS
# 2. Ones for Reid--simpler (pvals with sign)
# 3. Comparison of DMRseq and DSS
# 4. At regions of interest...

# CONSTANTS
GEN <- "hg38"
FS <- 16 # fontsize
N <- 5

#Colors for plotting
NCOL <- 7
colors <- pal_nejm("default")(NCOL)

# Parser
parser <- ArgumentParser()
parser$add_argument("--chr", default= "chr6", help='Chromosome to generate figures on')
parser$add_argument("--regions_file", default= "./regions-of-interest.csv", help='CSV specifying other regions to plot')
parser$add_argument("--array_dmps_file", default = "pub-2018Madrid-DMPs.lifted.csv")
parser$add_argument("--phenotypes_file", default = "../../data/meta/phenos-cleaned.csv")
parser$add_argument("--odir", default = "/figs/")
args <- parser$parse_args()

# Lazy way to recode variables
CHR <- args$chr
idir <- paste0("./result-", CHR, "/")
odir <- paste0(file.path(idir, "figs/"))
dir.create(odir, showWarnings=FALSE)

# More data
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) <- CHR
TxDb(Homo.sapiens) <- txdb
tx.hs <- unlist(transcriptsBy(Homo.sapiens, columns = "SYMBOL", by = "gene"))
ix <- seqnames(tx.hs) == CHR
tx.hs <- tx.hs[ix, ]


# Read in data
roi.df <- read_csv(args$regions_file, show_col_types = FALSE)
Meth <- fread(file.path(idir, "methylation.csv"))
load(file.path(idir, "models.RData"))

# Call DMRs--allow user to input parameters from terminal
dmrs <- callDMR(test.cohort, p.threshold=0.01, dis.merge = 2500, pct.sig = 0.5, minCG = 3)

#--> Split into two groups
load.samples <- filt.df %>% filter(cohort == "AD") %>% pull(sample) %>% as.character()
control.samples <- filt.df %>% filter(cohort == "CONTROL") %>% pull(sample) %>% as.character()

prepare_granges <- function(DT, chr=CHR){
    # Munges dataframe (changes column names)
    # to be amenable to GRanges() constructor
    DT.2 <- DT
    
    colnames(DT.2)[colnames(DT.2) == "pos"] <- "start"
    DT.2$end <- DT.2$start + 2
    DT.2$chr <- chr
    return(DT.2)
}

make_gene_track <- function(start, stop, genome=GEN, chr=CHR){
    #out <- GeneRegionTrack(txdb, chromosome = CHR, from = start, to = stop, fontsize = FS, 
    #                        name = "Reference Genes", fill = colors[7], transciptAnnotation = "symbol")

    biomTrack <- BiomartGeneRegionTrack(genome = genome, chromosome = chr, 
                    name = "ENSEMBL", start = start, end = stop, 
                    transciptAnnotation = "symbol", 
                    fontsize = FS, fill = colors[7], showID = TRUE)
    return(biomTrack)
}

Meth.2 <- prepare_granges(Meth)
group <- filt.df$cohort[match(colnames(Meth)[-1], filt.df$sample)]

# "Simple" tracks without too many parameters
dmrs.track <- AnnotationTrack(GRanges(dmrs), genome = GEN, name = "DMRs (DSS)", 
                            fontsize = FS, fill = colors[3], col = "white")
ideogram.track <- IdeogramTrack(genome = GEN, chromosome = CHR, fontsize = FS - 2)
axis.track <- GenomeAxisTrack()


tmp <- data.frame(start = test.cohort$pos, end = test.cohort$pos + 2, 
                chr = rep(CHR, nrow(test.cohort)), p = -1*log10(test.cohort$pvals))
p.track <- DataTrack(GRanges(tmp), type=c("p", "smooth"), baseline= -log10(0.01), col=colors[1],
                       col.baseline="black", name="-log10(p)", fontsize = FS)

# Tracks from custom functions
diff.track <- make_difference_track(Meth.2, load.samples, control.samples,
                                genome=GEN, chr=CHR, color=colors[6], fontsize=FS)
top.list <- list(ideogram.track) # ontop of genomic information
bottom.list <- list(dmrs.track, p.track, diff.track, axis.track) # below genomic information


create_output_name <- function(odir, nCG, chr, start, pad){
    # odir := output directory
    # nCG := number of CpGs in DMR
    # chr := chromosome 
    # start := starting position of DMR
    bn <- paste0(as.character(nCG), "-", chr, "-start-", as.character(start),"-pad-",as.character(pad), ".pdf")
    return(file.path(odir, bn))
}


plot_and_save <- function(start, stop, tracks.list, ofile){
    # start : starting position
    # stop: stopping position
    # tracks.list: list of tracks to plot, produced beforehand
    # ofile: output
    pdf(ofile)
    plotTracks(tracks.list, from = start, to = stop, background.title = colors[6], 
                collapseTranscripts = "meta", showID = TRUE, transcriptAnnotation = "symbol")
    dev.off()
}

get_gene_track <- function(start, stop, pad){
    a <- start - max(pad)
    b <- stop + max(pad)
    return(make_gene_track(a, b))
}

plot_by_ix <- function(ix, pad=c(1000, 5000, 1000), top=top.list, bottom=bottom.list){
    a <- dmrs$start[ix]
    b <- dmrs$end[ix]
    nCG <- dmrs$nCG[ix]

    gene.track <- get_gene_track(a,b, pad)
    all.tracks <- append(append(top.list, gene.track), bottom.list)

    for (p in pad){
        ofile <- create_output_name(odir, nCG, CHR, a, p)
        plot_and_save(a-p, b+p, all.tracks, ofile)
        #plot_and_save(a-p, b+p, gene.track, ofile)
    }
}


lapply(1:N, plot_by_ix)


#--> Regions of interest outside of DMRs
roi.df <-roi.df %>% filter(chr == CHR)
p <- 1000
nCG <- 0

for (ix in 1:nrow(roi.df))
{
    a <- roi.df$start[ix]
    b <- roi.df$end[ix]

    odir.tmp <- file.path(odir, roi.df$name[ix])
    dir.create(odir.tmp, showWarnings=FALSE)

    gene.track <- get_gene_track(a,b, p)
    all.tracks <- append(append(top.list, gene.track), bottom.list)
    ofile <- create_output_name(odir.tmp, nCG, CHR, a, p)
    plot_and_save(a-p, b+p, all.tracks, ofile)
}


format_locus <- function(seq, a, b){
    out <- paste0(seq, ":", as.character(a), "-", as.character(b))
}


window <- 10000
# Nearest gene
dmrs.gr <- GRanges(dmrs[1:N, ]) 
overlap <- findOverlaps(dmrs.gr + window, unstrand(tx.hs))

tx.hs.shrunk <- tx.hs[subjectHits(overlap), ]
dmrs.gr.expanded <- dmrs.gr[queryHits(overlap), ]

dist <- distanceToNearest(dmrs.gr.expanded, tx.hs.shrunk)
dd <- as.data.frame(dist)$distance

dmrs.df <- as.data.frame(dmrs.gr.expanded, row.names = NULL) %>%
            transmute(dmr.locus = format_locus(seqnames, start, end))

dmrs.df$distance <- dd

ann.df <- as.data.frame(tx.hs.shrunk[, "SYMBOL"], row.names = NULL) %>%
            transmute(gene.locus = format_locus(seqnames, start, end),
             gene.strand = strand, gene.name = SYMBOL)

comb.df <- cbind(dmrs.df, ann.df)

write_csv(comb.df, file = file.path(idir, "closest-genes.csv"))

# TEST
biomTrack <- BiomartGeneRegionTrack(genome = GEN, chromosome = CHR, 
                    name = "ENSEMBL", start = a, end = b,
                    fontsize = FS, fill = colors[7])


png("test.png")
plotTracks(biomTrack, collapseTranscripts = "meta", showID = TRUE, transcriptAnnotation = "symbol")
dev.off()