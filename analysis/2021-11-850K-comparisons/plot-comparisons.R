
library(tidyverse)
library(GenomicRanges)
library(liftOver)
library(wiscR)
library(parallel)
library(tools)
library(data.table)
library(viridis)

DATE <- (Sys.time() %>% str_split(" "))[[1]][1]

ROOT_DIR <- "../../data/"
enzyme_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].tsv", recursive = TRUE, full = TRUE)
array_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].bed", recursive = TRUE, full = TRUE)
chain <- import.chain("./hg19ToHg38.over.chain")

cobble_and_plot <- function(enzyme.file){
# enzyme.file is the relative path to the enzymatic data (whole-genome)
# array.file is the relative path to the EPIC 850K array file 
    base.name <- file_path_sans_ext(basename(enzyme.file))
    ix <- str_detect(array_files, base.name)

    if (all(!ix)) { 
        return(base.name) 
    } else {
        array.file <- array_files[ix]
    }

    # enzymatic (whole gneome data)
    enzyme.df <-fread(enzyme.file) %>% 
        dplyr::mutate(wg_methylation = methylated / (methylated + unmethylated)) %>% 
        dplyr::mutate(strand = ifelse(strand == "+-", "*", strand)) %>% 
        dplyr::mutate(stop = pos + 1) 
    enzyme.gr <- makeGRangesListFromDataFrame(enzyme.df, seqnames.field = "chr", start.field = "pos", keep.extra.columns = TRUE)

    array.df <- fread(array.file)
    
    # liftover chain can be downloaded thru the loaded package or via UCSC and command line
    # TODO link to this    
    # genomic ranges is where the magic happens
    array.gr <- GRanges(seqnames = array.df$chrom,
                        ranges = array.df$chromStart ,
                        array_methylation = array.df$Beta)
    # lifting 800k array sites to whole genome is faster than the other way around
    lifted.gr <- liftOver(array.gr, chain = chain)
    overlap.enzyme <- findOverlaps(enzyme.gr, lifted.gr, type = "start")

    df <- data.frame(lifted.gr[subjectHits(overlap.enzyme), 'array_methylation'],
                    enzyme.gr[queryHits(overlap.enzyme), 'wg_methylation']) %>% 
        dplyr::select(c("array_methylation", "wg_methylation")) %>% 
        drop_na()

    N <- nrow(df)
    rho <- cor(df$array_methylation, df$wg_methylation)

    p <- df %>% 
        ggplot(aes(x = array_methylation, wg_methylation)) +
        geom_hex(aes(fill=log(..count..))) +
        xlab("EPIC 850K Array Methylation") +
        ylab("NEB Next Methylation") +
        ggtitle(paste("Correlation for sample", base.name, "is", round(rho, 3), 
                "\nwith", round(N, -3), "overlapping loci")) +
        scale_fill_viridis(option = "magma") +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        geom_abline(intercept = 0, slope = 1)

    q <- p + wiscR::light_theme() 
    wiscR::save_plot(q, paste0("./figs/", base.name, "_corr.png"))
}

mclapply(enzyme_files, cobble_and_plot, mc.cores = 8)
