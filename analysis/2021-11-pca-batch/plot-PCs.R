library(tidyverse)
library(wiscR)
library(argparse)
library(viridis)

#--> Samplesheet with columns like 'sample', 'machine', 'pool'
samples.df <- read.table("../../data/meta/meta-data.tsv", header=TRUE) %>%
                dplyr::mutate(pool_group = paste(pool, group, sep="-"))

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('-idir', default = "../../data/prin_comps/", help = "input directory.")
parser$add_argument('-odir', default = paste0("figs-", Sys.Date()), help='where to store output figures.')
args <- parser$parse_args()
dir.create(args$odir, showWarnings = FALSE)

my_files <- list.files(path = args$idir, pattern =  ".tsv", recursive = TRUE, full = TRUE)

plot_pca <- function(df, color_by_str, label = FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
    p <- df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        ggtitle("First two PCs of methylation matrix") 
        #labs(caption = paste0("N observations: ", n, "\t", "N missing: ", n_na, "\n", "N imputed w/ local mean: ", n_imputed))

    if (length(unique(df[[color_by_str]])) > 10){
        p <- p + scale_color_viridis(option = "magma")
    }

    if (label){
        p <- p + geom_text(aes(label=sample),hjust=0, vjust=0, size = 16)
        prepend <- paste0(prepend, "_labeled")
    }

    ofile <- file.path(args$odir, paste0(prepend, "_", color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}


for (ff in my_files){
    prepend <- str_replace(basename(ff), ".tsv", "")
    df <- read_tsv(ff) %>% left_join(samples.df, by = "sample")

    plot_pca(df, "machine")
    plot_pca(df, "pool")
    plot_pca(df, "num_null")
    plot_pca(df, "num_not_null")
    plot_pca(df, "mean_methylation")
    plot_pca(df, "mean_methylation", label=TRUE)
    plot_pca(df, "mean_coverage")
    plot_pca(df, "median_coverage")
    plot_pca(df, "batch")
    print(ff)
}
