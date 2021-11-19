#TODO: get PCA working on five samples
#Go chromosome by chromosome? Or can we do all at once.
#
library(data.table)
library(tidyverse)
library(ggfortify)
library(gmodels)

#TODO: make it command line? is this worth the extra dependencies?
#TODO: take multiple "ROOT_DIR"s

ROOT_DIR <- "../../data/2021-11-03-batch01/"
samples.df <- read.table(paste0(ROOT_DIR, "meta-data.tsv"), header=TRUE)


for (row.ix in 1:nrow(samples.df)){
    tmp <- samples.df[row.ix, ]
    file_name <- file.path(ROOT_DIR, paste(tmp$pool,tmp$group, sep= "-"), "04-extract", paste0(tmp$sample, ".tsv"))
    if (file.exists(file_name)){
        print("hell yeah")
    }
}

read_and_filter <- function(file_name, chr_vec=c("chr1")){
    sample_name <- str_remove(file_name, ".tsv")

    DT <- fread(file_name)
    DT.filtered <- DT[intersect(which(chr %in% chr_vec), which(methylated + unmethylated > 0)), ]
    
    # Wont work with multiple chromosomes as is...
    DT.filtered$tmp <- DT.filtered$methylated / (DT.filtered$methylated + DT.filtered$unmethylated)
    DT.sub <- DT.filtered[, c("tmp", "pos")]
    names(DT.sub)[1] <- sample_name
    return(DT.sub)
    #Consider--if we only have negative strand, go ahead and cast it to postivie strand, move pos back one
}

# Do a full join?
my_files <- c("100.tsv", "101.tsv", "104.tsv", "105.tsv", "106.tsv")
is_first <- TRUE

for (ff in my_files){
    if (is_first){
        DT <- read_and_filter(ff)
        is_first <- FALSE
    } else {
        DT <- merge(DT, read_and_filter(ff), by = "pos")
    }
}

DTT <- transpose(DT[, -c("pos")])

pca.big <- fast.prcomp(DTT, center = TRUE, scale = FALSE)

p <- autoplot(pca.big)
ggsave("test-big.png", p)
