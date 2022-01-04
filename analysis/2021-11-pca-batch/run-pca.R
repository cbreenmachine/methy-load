#TODO: get PCA working on five samples
#Go chromosome by chromosome? Or can we do all at once.

library(data.table)
library(tidyverse)
library(gmodels)
library(wiscR)
library(argparse)

DATE <- (Sys.time() %>% str_split(" "))[[1]][1]

samples.df <- read.table("../../data/meta/meta-data.tsv", header=TRUE) %>%
                dplyr::mutate(pool_group = paste(pool, group, sep="-"))

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--pca', action='store_true', help='load raw files and recompute PCA')
parser$add_argument('--quick', action='store_true', help='run on one sec chrom only to see if everything works as expected')

args <- parser$parse_args()


ROOT_DIR <- "../../data/2021-12-02-methylation-coverage/"
my_files <- list.files(path = ROOT_DIR, pattern =  "meth_cov.tsv", recursive = TRUE, full = TRUE)

if (args$quick){ my_files <- my_files[length(my_files)]}

print(my_files)

for (ff in my_files){
    # Gets chr12 from string
    chr_name <- (basename(ff) %>% str_split("_"))[[1]][1]
    prepend <- paste0(DATE, "-", chr_name, "-")



    # READ data
    DT <- fread(ff)

    DT[ , M := methylated / coverage]
    DT.wide <- dcast(DT, sample~pos, value.var="M")
    DT.long <- dcast(melt(DT.wide, id.vars="sample", variable.name = "pos"), pos ~ sample)


    mask <- frollmean(DT.long[, -c("pos")], n = 20, fill = 0, na.rm = TRUE)
    setDT(mask)

    print(object.size(x=lapply(ls(), get)), units="Gb")

    for (j in seq_len(ncol(DT.wide))[-1] ){
        ix_range <- which(is.na(DT.wide[[j]]))
        set(DT.wide, ix_range, j, as.numeric(unlist(mask[j-1, ..ix_range])))
    }


    n <- nrow(DT.long) * (ncol(DT.long) - 1)
    n_na <- sum(sum(is.na(DT.long)))
    n_after <- sum(sum(is.na(DT.wide)))
    n_imputed <- n_na - n_after


    # Lost the positions here
    DT.wide <- nafill(DT.wide, fill = 0)
    setDT(DT.wide)
    names(DT.wide) <- c("sample", as.character(DT.long$pos))


    pca.big <- fast.prcomp(DT.wide[ , -c("sample")], center = TRUE, scale = FALSE) 

    df <- as.data.frame(pca.big$x)
    df$sample <- as.numeric(DT.wide$sample)
    df <- dplyr::inner_join(df, samples.df, by="sample")

    s <- pca.big$sdev # extract standard deviation
    v <- c(0, s^2) # zero-pad and square for variance
    percent_v <- round(diff( cumsum(v) / sum(v)) * 100, 2) # get % variance explained

    write.csv(df, paste0(prepend, "pca.csv"), row.names = FALSE, quote = FALSE)
    write.csv(percent_v, paste0(prepend, "percent_var_exp.csv"), row.names = FALSE, quote = FALSE)



    plot_pca <- function(df, color_by_str){
        # 2d scatterplot of first two principal components
        # colored by 'color_by_str' which should be string
        p <- df %>%
            ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
            geom_point(size = 6, alpha = 0.7) +
            wiscR::light_theme() +
            xlab(paste0("PC1 (", percent_v[1], "%)")) +
            ylab(paste0("PC2 (", percent_v[2], "%)")) +
            ggtitle("First two PCs of methylation matrix") +
            labs(caption = paste0("N observations: ", n, "\t", "N missing: ", n_na, "\n", "N imputed w/ local mean: ", n_imputed))

        wiscR::save_plot(p, paste0(prepend, color_by_str, ".png"))
    }

    
    plot_pca(df, "machine")
    plot_pca(df, "pool")
    plot_pca(df, "pool_group")

    rm(DT.wide, DT.long, mask)
    gc()

}