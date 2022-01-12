
library(tidyverse)
library(GenomicRanges)
library(liftOver)
library(wiscR)
library(parallel)
library(tools)
library(data.table)
library(viridis)

DATE <- Sys.Date()
odir <- paste0("./", DATE, "-figs/")
dir.create(odir)

ROOT_DIR <- "../../data/"
enzyme_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].tsv", recursive = TRUE, full = TRUE)
enzyme_samples <- str_remove(basename(enzyme_files), ".tsv")
array_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].bed", recursive = TRUE, full = TRUE)
array_samples <- str_remove(basename(array_files), ".bed")

valid <- intersect(array_samples, enzyme_samples)

enzyme_files <- enzyme_files[which(enzyme_samples %in% valid)]
array_files <- array_files[which(array_samples %in% valid)]

input.df <- data.frame(enzyme.file = enzyme_files, array.file = array_files, base.name = valid)
print(input.df)
chain <- import.chain("./hg19ToHg38.over.chain")

plot_hex_corr <- function(df, base.name){
    # makes geom_hex plot with array methylaytion on x axis and 
    # whole genome methylation on y.

    # Summary statistics first
    N <- nrow(df) # number of loci in both techs
    rho <- cor(df$array_methylation, df$wg_methylation) 
    
    p <- df %>% 
        ggplot(aes(x = array_methylation, wg_methylation)) +
        geom_hex(aes(fill=log(..count..))) +
        xlab("EPIC 850K Array Methylation") +
        ylab("NEB Next Methylation") +
        ggtitle(paste("Correlation for sample", base.name, "is", round(rho, 3), 
                "\nwith", round(N, -3), "overlapping loci")) +
        scale_fill_viridis(option = "magma", breaks = c(0,5,10)) +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        geom_abline(intercept = 0, slope = 1, size = 1.2) +
        wiscR::light_theme() 

    return(p)

}


cobble_and_plot <- function(ix){
# enzyme.file is the relative path to the enzymatic data (whole-genome)
# array.file is the relative path to the EPIC 850K array file 
    base.name <- input.df[ix, ]$base.name
    array.df <- fread(input.df[ix, ]$array.file)
    # enzymatic (whole gneome data)
    enzyme.df <-fread(input.df[ix, ]$enzyme.file) %>% 
        dplyr::mutate(wg_methylation = methylated / (methylated + unmethylated)) %>% 
        dplyr::mutate(strand = ifelse(strand == "+-", "*", strand)) %>% 
        dplyr::mutate(stop = pos + 1) 
    enzyme.gr <- makeGRangesListFromDataFrame(enzyme.df, seqnames.field = "chr", start.field = "pos", keep.extra.columns = TRUE)
    
    # liftover chain can be downloaded thru the loaded package or via UCSC and command line
    # TODO link to this    
    # genomic ranges is where the magic happens
    array.gr <- GRanges(seqnames = array.df$chrom,
                        ranges = array.df$chromStart ,
                        array_methylation = array.df$Beta)
    # lifting 800k array sites to whole genome is faster than the other way around
    lifted.gr <- liftOver(array.gr, chain = chain)
    overlap.enzyme <- findOverlaps(enzyme.gr, lifted.gr, type = "start")

    # Packs information into common dataframe
    df <- data.frame(lifted.gr[subjectHits(overlap.enzyme), 'array_methylation'],
                    enzyme.gr[queryHits(overlap.enzyme), 'wg_methylation']) %>% 
        dplyr::select(c("array_methylation", "wg_methylation")) %>% 
        drop_na()

    
    # Standsrd (full-range)
    p <- plot_hex_corr(df, base.name)
    wiscR::save_plot(p, paste0(odir, base.name, "_corr.png"))

    # Only the middle bits
    left <- 0.3
    right <- 0.7

    df.2 <- df %>% 
            dplyr::filter(left < array_methylation & array_methylation < right) %>%
            dplyr::filter(left < wg_methylation & wg_methylation < right)

    q <- plot_hex_corr(df.2, base.name) + 
        scale_fill_viridis(option = "magma") + # cannot manually specify breaks
        xlim(c(left-0.05,right+0.05)) +
        ylim(c(left-0.05,right+0.05))
    wiscR::save_plot(q, paste0(odir, base.name, "_mid_corr.png"))
       
}


mclapply(1:nrow(input.df), cobble_and_plot, mc.cores = 6)
