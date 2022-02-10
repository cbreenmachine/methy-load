
library(DSS)
library(Gviz)
library(tidyverse)
library(data.table)
library(ggsci)
library(parallel)
library(scales)
library(wiscR)
library(argparse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# user inputs
GEN <- "hg38"
FS <- 18 # fontsize


parser <- ArgumentParser()
parser$add_argument("--chr", default= "chr6", help='Chromosome to generate figures on')
args <- parser$parse_args()

# Lazy way to recode variables
CHR <- args$chr

# supp.df <- read_csv()
idir <- paste0("./result-", CHR, "/")

roi.df <- read_csv("regions-of-interest.csv")
M <- fread(file.path(idir, "M.csv"))
Cov <- fread(file.path(idir, "Cov.csv"))
Meth <- fread(file.path(idir, "methylation.csv"))
load(file.path(idir, "test-values.RData"))
array.df <- read_csv("pub-2018Madrid-DMPs.lifted.csv") %>%
            filter(seqnames == CHR)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#Colors
NCOL <- 7
colors <- pal_nejm("default")(NCOL)


df <- read_csv("../../data/meta/phenos-cleaned.csv")
valid_samples <- intersect(colnames(Meth)[-1], df$sample)
filt.df <- df %>% filter(sample %in% valid_samples)

filt.df %>%
    group_by(cohort) %>%
    summarize(count = n()) %>%
    write_csv(file.path(idir, "cohort-counts.csv"))

filt.df %>%
    group_by(sex) %>%
    summarize(count = n()) %>%
    write_csv(file.path(idir, "sex-counts.csv"))



# Call DMRs--allow user to input parameters from terminal
dmrs <- callDMR(test.cohort, p.threshold=0.01, dis.merge = 2500, pct.sig = 0.5, minCG = 3)
write_csv(dmrs, file.path(idir, "dmrs.csv"))


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


make_difference_track <- function(Meth, load.samples, control.samples, chr = CHR){
    
    # Compute mean difference
    meandiff <- rowMeans(Meth[ , load.samples, with = FALSE], na.rm = TRUE) - 
            rowMeans(Meth[ , control.samples, with = FALSE], na.rm = TRUE)
    # Construct GRanges
    dt <- data.table(diff=meandiff, chr=CHR, start=Meth$start, end=Meth$end) %>%
            drop_na() # important to prevent errors
    gr <- GRanges(dt)
    track <- DataTrack(gr, type=c("p", "smooth"), baseline=0, col=colors[2],
                       col.baseline="black", name="Mean Methylation Difference", 
                       fontsize = FS)
    return(track)
} 


make_methylation_track <- function(Meth, group){
    # Meth matrix should have been run through 
    DataTrack(GRanges(Meth), 
              groups = group, 
              name = "Methylation",                             
              type = c("a", "p"), 
              legend = TRUE) %>%
    return()
}


make_cpg_track <- function(start, stop){
    cpgIslands <- UcscTrack(genome = GEN, chromosome = CHR, 
                        track = "cpgIslandExt", 
                        from = start, to = stop,
                        trackType = "AnnotationTrack", 
                        start = "chromStart", end = "chromEnd", 
                        id = "name", shape = "box", fill = colors[4], 
                        name = "CpG Islands", fontsize = FS)
    return(cpgIslands)
}


make_snp_track <- function(start, stop){
    
    snpLocations <-  UcscTrack(genome = GEN, chromosome = CHR, 
                           track = "snp141", 
                           from = start, to = stop,
                           trackType = "AnnotationTrack", 
                           start = "chromStart", end = "chromEnd", 
                           id = "name", feature = "func", 
                           strand = "strand", shape = "box", 
                           stacking = "dense", fill = colors[2],
                           name = "SNPs", fontsize = FS)
    return(snpLocations)
}

make_p_track <- function(test.cohort){
     # Meth matrix should have been run through 

    data.frame(chr = test.cohort$chr, start = test.cohort$pos, stop = test.cohort$pos + 2, 
               log.p = -1*log10(test.cohort$pvals)) %>%
        drop_na() %>%
        GRanges() %>%
        DataTrack(name = "-log10(p)", type = c("smooth", "p"), 
              legend = TRUE, col=colors[1], fontsize = FS) %>%
    return()
}


make_array_p_track <- function(array.df){
    GRanges(array.df) %>%
    DataTrack(name = "-log10(p) (Array)", type = "p", 
              legend = TRUE, col=colors[2], fontsize = FS) %>%
    return()
}

make_gene_track <- function(start, stop, chr=CHR){
    out <- GeneRegionTrack(txdb, chromosome = CHR, from = start, to = stop, fontsize = FS, name = "Ref Genes",
                           transcriptAnnotation = "symbol", fill = colors[7])
    return(out)
} 

Meth.2 <- prepare_granges(Meth)
group <- filt.df$cohort[match(colnames(Meth)[-1], filt.df$sample)]


# "Simple" tracks without too many parameters
dmrs.track <- AnnotationTrack(GRanges(dmrs), genome = GEN, name = "DMRs", fontsize = FS,
                             fill = colors[3], col = "white")
ideogram.track <- IdeogramTrack(genome = GEN, chromosome = CHR, fontsize = FS - 2)
#axis.track <- GenomeAxisTrack()
p.track <- make_p_track(test.cohort)
#array.p.track <- make_array_p_track(array.df)

# Tracks from custom functions
diff.track <- make_difference_track(Meth.2, load.samples, control.samples)
# Show up on top of genomic information
top.list <- list(ideogram.track)
# Show up below gneomic info
bottom.list <- list(dmrs.track, p.track, diff.track)



plot_all_tracks <- function(start, stop, top.tracks.list, bottom.tracks.list, 
                            odir, pad=c(1000, 5000, 10000)){
    # start : starting position
    # stop: stopping position
    # tracks_list: list of tracks to plot, produced beforehand
    # odir: where to save pngs
    # pad: vector specifying how big a window
    
    # File setup
    dir.create(odir, showWarnings=FALSE)
    
    # Make the genes track based on the largest window
    a <- start - max(pad)
    b <- stop + max(pad)
    
    # All the tracks
    ref.genes.track <- make_gene_track(start = a, stop = b)
    #snp.track <- make_snp_track(start = a, stop = b)
    cpg.track <- make_cpg_track(start = a, stop = b)
    
    # Cobble together
    mid.tracks.list <- list(ref.genes.track, cpg.track)
        
    # Loop thru padding values
    for (p in pad){
        a <- start - p
        b <- stop  + p
        
        tracks.list <- append(top.tracks.list, mid.tracks.list)
        tracks.list <- append(tracks.list, bottom.tracks.list)
        
        
        # Save, how do we feel about default sizes?
        pdf(paste0(odir, as.character(start), "-", as.character(p), "-DMRs.pdf"))
        plotTracks(tracks.list, from = a, to = b, 
                   background.title = colors[6])
        dev.off()
    }    
}



plot_by_ix <- function(ix){
    plot_all_tracks(dmrs$start[ix], dmrs$end[ix], top.list, bottom.list, file.path(idir, "/figs/"))
}

mclapply(1:nrow(dmrs), plot_by_ix, mc.cores = 4)


#--> Regions of interest outside of DMRs
roi.df <-roi.df %>% filter(chr == CHR)


for (ix in 1:nrow(roi.df)){
    plot_all_tracks(roi.df$start[ix], roi.df$end[ix], top.list, bottom.list, file.path(idir, "figs/", roi.df$name[ix], "/"), pad = c(1000, 2000))
}
