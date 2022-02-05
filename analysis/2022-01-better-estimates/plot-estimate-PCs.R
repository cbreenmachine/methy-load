library(tidyverse)
library(wiscR)
library(argparse)
library(viridis)

# write covariates to a file in odir
parser <- ArgumentParser(description='Methylation estimates')
parser$add_argument('--ifile', default ="./data/PC-cell-compositions.csv", help='where to store output data.')
parser$add_argument('--odir', default ="./figs/PC-cell-compositions/")
args <- parser$parse_args()

dir.create(args$odir, showWarnings = FALSE)

#--> Samplesheet with columns like 'sample', 'machine', 'pool'
samples.df <- read_csv("../../data/meta/phenos-cleaned.csv", col_types = cols()) %>% 
                dplyr::mutate(pool_group = paste(pool, group, sep="-"))
other.pc.df <- read_tsv("../../data/prin_comps/PC_chr22.tsv", col_types = cols()) %>%
                dplyr::select(-starts_with('PC'))
samples.df <- merge(samples.df, other.pc.df, by = "sample")


plot_pca <- function(df, color_by_str, label = FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
    #TODO: clean label string (title case and _ to " ")
    #TODO: 3 breaks in legend--min, mean, max? Why are they inconsistent
    p <- df %>%
        drop_na() %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 6, alpha = 0.8) +
        wiscR::light_theme() +
        ggtitle("First two PCs of methylation matrix") 

    if ((length(unique(df[[color_by_str]])) > 10) & (is.numeric(df[[color_by_str]]))) {
        p <- p + scale_color_viridis(option = "magma")
    }

    if (label){
        # adds '101', '432', etc. to the points
        # font sizes and colors are a bit wacky, but helps identify outliers
        p <- p + geom_text(aes(label=sample),hjust=0.5, vjust=0.1, size = 16) 
        color_by_str <- paste0(color_by_str, "_labeled")
    }

    ofile <- file.path(args$odir, paste0(color_by_str, ".png"))
    wiscR::save_plot(p, ofile)
}

df <- read_csv(args$ifile, col_types = cols() ) %>% 
        merge(samples.df, by = "sample")

all_vars <- names(df)
good_vars <- all_vars[!str_detect(all_vars, "PC|sample")]

for (vv in good_vars){
    plot_pca(df, vv)
    print(vv)
}

plot_pca(df, "mean_methylation", label=TRUE)

