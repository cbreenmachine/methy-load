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
    library(biomaRt)
    library(methanalyzeR)
}) 



# CONSTANTS
ALPHA <- 0.01
GEN <- "hg38"
FS <- 16 # fontsize
N <- 5

#Colors for plotting
NCOL <- 7
colors <- pal_nejm("default")(NCOL)

# Parser
parser <- ArgumentParser()
parser$add_argument("--chr", default= "chr6", help='Chromosome to generate figures on')
parser$add_argument("--experiment", default= "/smooth-150-PCs-2/", help='Chromosome to generate figures on')
parser$add_argument("--regions_file", default= "./regions-of-interest.csv", help='CSV specifying other regions to plot')
parser$add_argument("--phenotypes_file", default = "../../data/meta/phenos-cleaned.csv")
parser$add_argument("--smoothing_window", default = 150)
parser$add_argument("--odir", default = "/figs/")
args <- parser$parse_args()


print("Parsed args")

# Lazy way to recode variables
CHR <- args$chr
idir <- file.path(paste0("./", CHR, "/"), args$experiment)
odir <- paste0(file.path(idir, "figs/"))
dir.create(odir, showWarnings=FALSE)

# Use the same biomart that we used for annotations
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Read in data
roi.df <- read_csv(args$regions_file, show_col_types = FALSE)
Meth <- fread(file.path(dirname(idir), "methylation.csv")) # %>% drop_na()
load(file.path(idir, "models.RData"))

# Smooth 
vv <- names(Meth)[-1] 
DF <- as.data.frame(Meth)
out <- lapply(vv, function(x) frollmean(DF[ , x], n = 150, fill = 0.8, align = "center", na.rm = T))
Meth.smooth <- data.table(do.call(cbind, out))
colnames(Meth.smooth) <- vv
Meth.smooth$pos <- Meth$pos

dmrs <- read_csv(file.path(idir, "called-dmrs.csv"),show_col_types = FALSE) %>% 
        separate(col = dmrLocus, into = c("chr", "start", "end")) %>%
        mutate(start = as.numeric(start), end = as.numeric(end))
closet.genes <- read_csv(file.path(idir, "closest-genes.csv"),show_col_types = FALSE)
load.samples <- filt.df %>% filter(cohort == "AD") %>% pull(sample) %>% as.character()
control.samples <- filt.df %>% filter(cohort == "CONTROL") %>% pull(sample) %>% as.character()

print("Loaded and processed data")


#--> How many are postivie/negative
diffs <- rowMeans(Meth.smooth[,..load.samples], na.rm = T) - rowMeans(Meth.smooth[,..control.samples], na.rm = T)
diffs.df <- data.frame(diff = diffs, pos = Meth.smooth$pos) 

.tally_signs <- function(df = diffs.df, start, end){
    signs <- df %>% filter(pos >= start) %>% filter(pos <= end) %>% pull(pos) %>% sign()
    nHyperMeth <- sum(signs == 1)
    return(nHyperMeth)
}

dmrs$nHyperMeth <- -1
#TODO: Vectorize
for (i in 1:nrow(dmrs)){
    dmrs$nHyperMeth[i] <- .tally_signs(diffs.df, dmrs$start[i], dmrs$end[i])
}


print("Tallied hyper methylated")
# tmp <- dmrs %>%
#     filter(str_detect(geneNames, "KIF25"))
# start <- tmp$start
# end <-  tmp$end

# Meth.smooth %>%
#     filter(pos >= start) %>% filter(pos <= end)


plot_dmr <- function(df=diffs.df, start, end, ofile){
    p <- df %>% 
        dplyr::filter(pos >= start) %>% 
        dplyr::filter(pos <= end) %>%
        ggplot(aes(x = pos, y = (diff))) +
        geom_point(size = 4, color = "orange") +
        geom_segment(aes(x=pos, xend=pos, y=0, yend=diff), size = 1.5, color = "grey") +
        wiscR::light_theme() +
        xlab("Genomic position") +
        ylab("Difference in mean methylation\n(LOAD - Cont)")
    wiscR::save_plot(p, ofile)
}

plot_all_dmrs <- function(dmrs, diffs.df){
    for (i in 1:nrow(dmrs)){
        start <- dmrs$start[i]
        end <- dmrs$end[i]

        ofile <- file.path(odir, paste0("dmr-", str_replace_all(dmrs$geneNames[i], "; ", "-"), ".png"))
        plot_dmr(diffs.df, start,end,ofile)
    }
}
plot_all_dmrs(dmrs, diffs.df)

print("Plotted all DMRs")

# HELPER FUNCTION ------------------------------------------------------------------------------
prepare_granges <- function(DT, chr=CHR){
    # Munges dataframe (changes column names)
    # to be amenable to GRanges() constructor
    DT.2 <- DT
    
    colnames(DT.2)[colnames(DT.2) == "pos"] <- "start"
    DT.2$end <- DT.2$start + 2
    DT.2$chr <- chr
    return(DT.2)
}

make_gene_track <- function(start, stop, genome=GEN, chr=CHR, biomart=ensembl){
    # Same data source as annotation
    biomTrack <- BiomartGeneRegionTrack(genome = genome, chromosome = chr, 
                    name = "ENSEMBL", start = start, end = stop, 
                    transciptAnnotation = "symbol", biomart = biomart,
                    fontsize = FS, fill = colors[7], showID = TRUE)
    return(biomTrack)
}


# Make common tracks ---------------------------------------------------------------------
Meth.2 <- prepare_granges(Meth.smooth)
group <- filt.df$cohort[match(colnames(Meth.smooth), filt.df$sample)]
group <- group[!is.na(group)]

# "Simple" tracks without too many parameters
dmrs.track <- AnnotationTrack(GRanges(dmrs), genome = GEN, name = "DMRs", 
                            fontsize = FS, fill = colors[3], col = "white")
ideogram.track <- IdeogramTrack(genome = GEN, chromosome = CHR, fontsize = FS - 2)
axis.track <- GenomeAxisTrack()


tmp <- data.frame(start = test.cohort$pos, end = test.cohort$pos + 2, 
                chr = rep(CHR, nrow(test.cohort)), p = -1*log10(test.cohort$pvals))
p.track <- DataTrack(GRanges(tmp), type=c("p", "smooth"), baseline= -log10(0.01), col=colors[1],
                       col.baseline="black", name="-log10(p)", fontsize = FS, genome = GEN)

# Tracks from custom functions
diff.track <- make_difference_track(Meth.2, load.samples, control.samples,
                                genome=GEN, chr=CHR, color=colors[6], fontsize=FS)
top.list <- list(ideogram.track) # ontop of genomic information
bottom.list <- list(dmrs.track, p.track, diff.track, axis.track) # below genomic information


create_output_name <- function(odir, nCG, chr, start, pad, geneNames){
    # odir := output directory
    # nCG := number of CpGs in DMR
    # chr := chromosome 
    # start := starting position of DMR
    bn <- paste0(as.character(nCG), "-", chr, "-start-", as.character(start),
                "-pad-", as.character(pad), "-", geneNames, ".pdf")
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

plot_by_ix <- function(ix, pad=c(1000, 5000, 10000), top=top.list, bottom=bottom.list){
    # For lapply, wrapper for make_gene_track and save...
    a <- dmrs$start[ix]
    b <- dmrs$end[ix]
    nCG <- dmrs$nCG[ix]
    geneNames <- str_replace_all(dmrs$geneNames[ix], "; ", "-")

    gene.track <- get_gene_track(a,b, pad)
    all.tracks <- append(append(top.list, gene.track), bottom.list)

    for (p in pad){
        ofile <- create_output_name(odir, nCG, CHR, a, p, geneNames)
        plot_and_save(a-p, b+p, all.tracks, ofile)
    }
}


lapply(1:nrow(dmrs), plot_by_ix)

print("Plotted all DMR tracks")

#--> Regions of interest outside of DMRs
roi.df <-roi.df %>% filter(chr == CHR)
p <- 1000
nCG <- 0

plot_roi_ix <- function(ix){
    a <- roi.df$start[ix]
    b <- roi.df$end[ix]

    # odir.tmp <- file.path(odir, roi.df$name[ix])
    # dir.create(odir.tmp, showWarnings=FALSE)

    gene.track <- get_gene_track(a,b, p)
    all.tracks <- append(append(top.list, gene.track), bottom.list)
    ofile <- create_output_name(odir, nCG, CHR, a, p, roi.df$name)
    plot_and_save(a-p, b+p, all.tracks, ofile)
}

lapply(1:nrow(roi.df), plot_roi_ix)


write_csv(dmrs, file.path(idir, "called-dmrs-w-hyper.csv")) 