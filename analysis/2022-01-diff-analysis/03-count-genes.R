suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(wiscR)
    library(argparse)
}) 


parser <- ArgumentParser()
parser$add_argument("--idir", default= "chr6/", help='Path to input models file')
parser$add_argument("--num_pcs", default= "2", help='Path to input models file')
parser$add_argument("--smoothing_ws", default= "100,150,200", help='Path to input models file')
parser$add_argument("--regions_file", default= "./regions-of-interest.csv", help='CSV specifying other regions to plot')
args <- parser$parse_args()

# TODO: add this function to some R package
valid.pcs <- as.vector(unlist(str_split(args$num_pcs, ",")))
valid.ws <- as.vector(unlist(str_split(args$smoothing_ws, ",")))

# Load data --------------------------------------------------------------------

all.files <- list.files(args$idir, pattern = "closest-genes.csv", recursive = TRUE, full = TRUE)
print(all.files)

load_data <- function(path){
    experiment <- basename(dirname(path))
    experiment.broken <- unlist(str_split(experiment, "-"))
    smoothing.ws <- experiment.broken[2]
    num.pcs <- experiment.broken[4]

    out.df <- read_csv(path,  col_types = cols()) %>%
                mutate(num.pcs = num.pcs, smoothing.ws = smoothing.ws, chr = basename(args$idir)) %>%
                select(c(chr, gene.name, num.pcs, smoothing.ws)) %>%
                unique()
    
    return(out.df)
}

df <- do.call("rbind", lapply(all.files, load_data)) %>%
        filter(num.pcs %in% valid.pcs) %>% filter(smoothing.ws %in% valid.ws)


count.df <- df %>% 
    unite(Params, c(num.pcs, smoothing.ws), sep = ":") %>%
    group_by(Params, gene.name) %>% 
    summarize(Count = n())
    

p <- count.df %>%
    ggplot(aes(x = Params, y = gene.name, fill = Count)) +
    geom_tile() +
    wiscR::light_theme() +
    xlab("Parameters (num PCs : smoothing window)") +
    ylab("Gene name") +
    scale_fill_gradient(low = "white", high = "#C5050C", limits=c(0,1)) +
    theme(legend.position="none") +
    labs(caption = "Red indicates that the given gene is within 10kb of a DMR.") +
    theme(axis.text.x = element_text(angle = 30))


ofile <- file.path(args$idir, "closest-genes-by-design.png")
wiscR::save_plot(p, ofile)