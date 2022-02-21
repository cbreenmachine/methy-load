suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(wiscR)
    library(argparse)
}) 


parser <- ArgumentParser()
parser$add_argument("--ifile", default= "chr6/smooth-250-PCs-0/models.RData", help='Path to input models file')
args <- parser$parse_args()


# Helper functions ---------------------------------------------------------------------
chain.file <- "hg19ToHg38.over.chain"
.check_chain <- function(){
    if (! file.exists(chain.file)){
        tmp <- file.path("https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/", chain.file)
        utils::download.file(tmp, destfile = paste0(chain.file, ".gz"))
        system(paste0("gzip -d ", chain.file, ".gz"))
    } 
}

lift_to_hg38 <- function(df){
    .check_chain()
    chain <- rtracklayer::import.chain(chain.file)

    # "Normal" coordinates, but renamed
    df <- df %>% 
            dplyr::rename(seqnames = chr) %>%
            dplyr::mutate(start = pos - 1, end = pos + 1)
    gr <- GenomicRanges::GRanges(df)
    gr.lifted <- rtracklayer::liftOver(gr, chain = chain)

    out.df <- as.data.frame(gr.lifted) %>% 
        dplyr::rename(chr = seqnames) %>% 
        dplyr::select(-c(group, group_name)) %>%
        dplyr::mutate(pos = start + 1) # put the pos back to sensible one-based
    return(out.df)
}


# Load data functions ---------------------------------------------------------------------

load(args$ifile)
DMPs.whole <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos, 
                         p.whole = test.cohort$pvals, adj.p.whole = -log10(test.cohort$pvals))

DMPs.array <- fread("../../../adrc-array/data/DMPs.csv") %>%
            filter(chr == DMPs.whole$chr[1]) %>%
            mutate(adj.p.array = -log10(P.Value), p.array = P.Value) %>%
            lift_to_hg38()

# Not neccessary since lifting figured out
print("Num positions in intersection")
length(intersect(DMPs.array$pos, DMPs.whole$pos))
print("Num positions in array / whole")
length(setdiff(DMPs.array$pos, DMPs.whole$pos))

# Cobble together by chromosome and position
join.df <- left_join(DMPs.whole, DMPs.array, by = c("chr", "pos"))
# Creaate variable telling us whether a location is in array loci or not
join.df <- join.df %>% 
            mutate(ArrayLocus = ifelse(is.na(adj.p.array), "NotInArray", "InArray")) %>%
            dplyr::select(-adj.P.Val)
roi.df <- read_csv("regions-of-interest.csv") %>% filter(chr %in% join.df$chr)

# Plot
plot_me <- function(joined.df, roi.df, ix=1, buffer=1000){
    p <- join.df %>%
        filter(pos > roi.df$start[ix] - buffer & pos < roi.df$end[ix] + buffer) %>%
        pivot_longer(cols = starts_with("adj.p"), names_to = "Tech", values_to = "adj.p") %>%
        ggplot(aes(x = pos, y = adj.p, color = ArrayLocus)) +
        geom_point(size = 4, alpha = 0.8) +
        xlab("Genomic position") +
        ylab("-log10(p) from DSS / whole-genome") +
        ggtitle(roi.df$name[ix]) +
        labs(caption = paste0("Range is +/-", as.character(buffer), " nt around ", roi.df$name[ix])) +
        wiscR::light_theme()
    odir <- (dirname(args$ifile))
    ofile <- file.path(odir, paste0(roi.df$name[ix], "-comp.png"))
    
    wiscR::save_plot(p, ofile)
}

plot_me(joined.df, roi.df)


rho <- cor(join.df$p.whole, join.df$p.array, use = "complete.obs")

p <- join.df %>%
    ggplot(aes(x = p.whole, p.array)) +
    geom_hex() +
    geom_rug(alpha = 0.1) +
    scale_fill_viridis(option = "magma", breaks = c(0,5,10)) +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    wiscR::light_theme() +
    ggtitle(paste0("Correlation: ", as.character(rho)))

ofile <- file.path(dirname(args$ifile), "pvalues-hex.png")
wiscR::save_plot(p, ofile)


p <- join.df %>%
    ggplot(aes(x = p.whole)) +
    geom_histogram() + 
    wiscR::light_theme() 
ofile <- file.path(dirname(args$ifile), "p-whole-genome.png")
wiscR::save_plot(p, ofile)



p <- join.df %>%
    ggplot(aes(x = p.array)) +
    geom_histogram() + 
    wiscR::light_theme() 
ofile <- file.path(dirname(args$ifile), "p-array.png")
wiscR::save_plot(p, ofile)