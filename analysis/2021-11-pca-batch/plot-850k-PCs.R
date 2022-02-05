library(data.table)
library(tidyverse)
library(wiscR)
library(parallel)
library(gmodels)
library(viridis)
# 153 deleted becaue it was missing lines...

idir <- "../../data/850k-array/"
odir <- "./figs-array"
prepend <- file.path(odir, "array_")
dir.create(odir, showWarnings = FALSE)

files <- list.files(idir, full = TRUE)


read_data <- function(file.name){
    DT <- fread(file.name, select = c("Beta"))
    names(DT)[1] <-str_remove(basename(file.name), ".bed")
    return(DT)
}

DT <- do.call(cbind, mclapply(files, read_data, mc.cores = 8))
DT.t <- (t(DT))

pca.big <- fast.prcomp(DT.t, center = TRUE, scale = FALSE)


tmp.1 <- read_csv("../../data/meta/phenos-cleaned.csv") %>% mutate(sample = as.character(sample))
tmp.2 <- read_tsv("../../data/prin_comps/PC_chr1.tsv") %>% select(!contains("PC")) %>% mutate(sample = as.character(sample))
tmp <- inner_join(tmp.1, tmp.2, by = "sample")
df <- as.data.frame(pca.big$x) %>% rownames_to_column("sample") %>%
        full_join(tmp, by = "sample")



plot_pca <- function(df, color_by_str, label=FALSE){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
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
        p <- p + geom_text(aes(label=sample),hjust=0, vjust=0, size = 12)
        prepend <- paste0(prepend, "labeled_")
    }
    
    wiscR::save_plot(p, paste0(prepend, color_by_str, ".png"))
}


for (vv in names(df)){
    print(vv)
    if (!str_detect(vv, "PC|sample")){
        plot_pca(df, vv)
    }
}

plot_pca(df, "mean_methylation", TRUE)



prepend <- file.path(odir, "array_no_outlier_")

df <- df %>% filter(sample != "101")

for (vv in names(df)){
    print(vv)
    if (!str_detect(vv, "PC|sample")){
        plot_pca(df, vv)
    }
}

plot_pca(df, "mean_methylation", TRUE)

