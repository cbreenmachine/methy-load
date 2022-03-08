suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(wiscR)
    library(argparse)
}) 


parser <- ArgumentParser()
parser$add_argument("--idir", default= "chr6/", help='Path to input models file')
parser$add_argument("--num_pcs", default= "0,2,4", help='Path to input models file')
parser$add_argument("--smoothing_ws", default= "0,50,100,250,500", help='Path to input models file')
parser$add_argument("--regions_file", default= "./regions-of-interest.csv", help='CSV specifying other regions to plot')
args <- parser$parse_args()

roi.df <- read_csv(args$regions_file, col_types = cols() ) %>% filter(chr== basename(args$idir))
# TODO: add this function to some R package
valid.pcs <- as.vector(unlist(str_split(args$num_pcs, ",")))
valid.ws <- as.vector(unlist(str_split(args$smoothing_ws, ",")))

# Load data --------------------------------------------------------------------

all.files <- list.files(args$idir, pattern = "models.RData", recursive = TRUE, full = TRUE)


load_data <- function(path){
    experiment <- basename(dirname(path))
    experiment.broken <- unlist(str_split(experiment, "-"))
    smoothing.ws <- experiment.broken[2]
    num.pcs <- experiment.broken[4]

    load(path)
    out.df <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos,
                        p = test.cohort$pvals, smoothing.ws = smoothing.ws, num.pcs = num.pcs)
    return(out.df)
}

df <- do.call("rbind", lapply(all.files, load_data)) %>%
        filter(num.pcs %in% valid.pcs) %>% filter(smoothing.ws %in% valid.ws)


ix <- 1
buffer <- 5000


p <- df %>%
    filter(roi.df$start[ix] - buffer < pos) %>% 
    filter(pos < roi.df$end[ix] + buffer) %>%
    mutate(num.pcs.2 = paste0(num.pcs, " PCs"), smoothing.ws.2 = paste0("Smooth ", smoothing.ws)) %>%
    ggplot(aes(x = pos, y = -log10(p))) +
    geom_point(size = 2) +
    geom_smooth(se=FALSE) + 
    facet_grid(vars(num.pcs.2), vars(smoothing.ws.2)) +
    wiscR::light_theme() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
    ggtitle(paste0("+/- ", as.character(buffer), " nt around ", roi.df$name[ix])) +
    xlab("Genomic position")
    
ofile <- file.path(args$idir, paste0(roi.df$name[ix], "-experiment-summary", ".png"))
wiscR::save_plot(p, ofile)