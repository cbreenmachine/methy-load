#TODO: get PCA working on five samples
#Go chromosome by chromosome? Or can we do all at once.
#

library(data.table)
library(tidyverse)
library(ggfortify)
library(gmodels)
library(wiscR)

#TODO: make it command line? is this worth the extra dependencies?
#TODO: take multiple "ROOT_DIR"s

ROOT_DIR <- "../../data/2021-11-03-batch01/"
my_files <- list.files(path = ROOT_DIR, pattern =  "[0-9][0-9][0-9].tsv", recursive = TRUE)

samples.df <- read.table(paste0(ROOT_DIR, "meta-data.tsv"), header=TRUE)


# TODO
#for (row.ix in 1:nrow(samples.df)){
#    tmp <- samples.df[row.ix, ]
#   file_name <- file.path(ROOT_DIR, paste(tmp$pool,tmp$group, sep= "-"), "04-extract", paste0(tmp$sample, ".tsv"))
#    if (file.exists(file_name)){
#        print("hell yeah")
#    }
#}

read_and_filter <- function(file_name, chr_vec=c("chr1")){
    sample_name <- str_remove(file_name, ".tsv") %>% basename()

    DT <- fread(file_name)
    DT.filtered <- DT[intersect(which(chr %in% chr_vec), which(methylated + unmethylated > 0)), ]
    
    # Wont work with multiple chromosomes as is...
    DT.filtered$tmp <- DT.filtered$methylated / (DT.filtered$methylated + DT.filtered$unmethylated)
    DT.sub <- DT.filtered[, c("tmp", "pos")]
    names(DT.sub)[1] <- sample_name
    return(DT.sub)
    #Consider--if we only have negative strand, go ahead and cast it to postivie strand, move pos back one
}


#included <- names(DT)[-1]
#subset.df <- samples.df %>%
#    dplyr::filter(sample %in% included) %>%
#    dplyr::arrange(-sample) %>%
#   dplyr::mutate(dummy = as.factor(1:n()))

# Order by number first when this automated more!!

is_first <- TRUE

for (ff in file.path(ROOT_DIR, my_files)){
    print(ff)
    if (is_first){
        DT <- read_and_filter(ff)
        print(paste("There are", nrow(DT), "in first file"))
        is_first <- FALSE
    } else {
        DT <- merge(DT, read_and_filter(ff), by = "pos", all=TRUE)
        print(paste("There are", nrow(DT), "after merging."))
    }
}

sample_means <- as.double(DT[, lapply(.SD, mean, na.rm = TRUE)][1,-1]) # drop position

DTT <- data.table::transpose(DT[, -c("pos")])
colnames(DTT) <- as.character(DT$pos)
rownames(DTT) <- names(DT)[-1]

for (j in seq_len(ncol(DTT))){ 
    ix <- which(is.na(DTT[[j]]))
    set(DTT, ix, j, sample_means[ix]) 
}

pca.big <- fast.prcomp(DTT, center = TRUE, scale = FALSE) 


df <- as.data.frame(pca.big$x)
df$sample <- as.numeric(rownames(DTT))
df <- dplyr::inner_join(df, samples.df, by="sample")

s <- pca.big$sdev # extract standard deviation
v <- c(0, s^2) # zero-pad and square for variance
percent_v <- round(diff( cumsum(v) / sum(v)) * 100, 2) # get % variance explained

plot_pca <- function(df, color_by_str){
    # 2d scatterplot of first two principal components
    # colored by 'color_by_str' which should be string
    p <- df %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color = color_by_str)) +
        geom_point(size = 5) +
        wiscR::light_theme() +
        xlab(paste0("PC1 (", percent_v[1], "%)")) +
        ylab(paste0("PC2 (", percent_v[2], "%)")) +
        ggtitle("First two PCs of methylation matrix")

    wiscR::save_plot(p, paste0(color_by_str, ".png"))
}
 
plot_pca(df, "machine")
plot_pca(df, "pool")
plot_pca(df, "group")

write.csv(df, "pca.csv", row.names = FALSE, quote = FALSE)
write.csv(percent_v, "percent_var_exp.csv", row.names = FALSE, quote = FALSE)

# May need to manually pull out matrix, and plot so that colors can be controlled
#ap <- autoplot(pca.big) + wiscR::light_theme()
#wiscR::save_plot(ap, "autoplot.png")
