
# System data, make sankey diagram
# Make a python crawler that creates the extracted file if it's not there already
suppressPackageStartupMessages( library(tidyverse) )

# max_possible <- 100 

# root_dir <- "../data/"
# raw_fastqs <- list.files(root_dir, pattern = "fastq.gz", recursive = TRUE) 
# raw_fastqs <- file.path(root_dir, raw_fastqs[!str_detect(raw_fastqs, "txt")])


# fq.df <-data.frame(file_name = basename(raw_fastqs),
#                     size = round(file.info(raw_fastqs)$size / 1024 ^ 3, 1)) %>%
#                     dplyr::mutate(sample = substr(file_name, 1, 3),
#                                   pg = str_extract(raw_fastqs, "pool[0-9][0-9]-group[0-9][0-9]"))

# fq.df %>%
#     group_by(sample) %>%
#     summarize(N = n(), min_size = min(size), max_size = max(size), pg) %>%
#     dplyr::arrange(N, min_size) %>% unique() %>% write_csv("file-sizes.csv")

# fq.df %>%
#     group_by(pg) %>%
#     summarize(N = n() / 2) %>% drop_na() %>% arrange(N) %>% write_csv("pool-group-sizes.csv")

# trimmed_fastqs <- list.files(root_dir, pattern = ".fq", recursive = TRUE)
# mapped <- list.files(root_dir, pattern = ".bam.csi", recursive = TRUE)
# called <- list.files(root_dir, pattern = "^[0-9][0-9][0-9].bcf.csi", recursive = TRUE)
# extracted <- list.files(root_dir, pattern = "^[0-9][0-9][0-9].tsv", recursive = TRUE) 

# N <- min(c(max_possible, length(raw_fastqs)))

# df <- data.frame(NumTrimmed = length(trimmed_fastqs) / 2,
#                  NumMapped  = length(mapped),
#                  NumCalled = length(called),
#                  NumExtracted = length(extracted))
# print(df)


count_raw <- function(path){
    return(length(list.files(file.path(path, "00-fastq"))) / 2)
}


count_trimmed <- function(path){
    tmp <- list.files(file.path(path, "01-fastq-trimmed"), pattern = "val_1.fq")
    return(length(tmp))
}

count_mapped <- function(path){
    tmp <- list.files(file.path(path, "02-mapping"), pattern = ".csi")
    return(length(tmp))
}

count_called <- function(path){
    tmp <- list.files(file.path(path, "03-calls"), pattern = ".md5")
    return(length(tmp))
}

dirs <- system('find ../data/ -maxdepth 2 -mindepth 2', intern = TRUE)
ix <- unlist(lapply(dirs, function(x) str_detect(x, pattern = "pool")))
dirs <- dirs[ix]


df <- data.frame(path = dirs) %>%
        mutate(poolgroup = basename(path), 
            nStart = unlist(lapply(path, count_raw)), 
            nTrimmed = unlist(lapply(path, count_trimmed)),
            nMapped = unlist(lapply(path, count_mapped)),
            nCalled = unlist(lapply(path, count_called))) %>%
        select(-path) %>%
        arrange(poolgroup)


ofile <- paste0(Sys.Date(), "-done.csv")
write_csv(df, file = ofile)


df <- read_tsv("../data/meta/ADRC-samplesheet-cleaned.tsv", show_col_types = FALSE)
extracted <- list.files("../data/", pattern = "^[0-9][0-9][0-9].tsv", recursive = TRUE) 
todo <- as.character(sort(setdiff(df$ALISCH_ID, str_remove(basename(extracted), ".tsv"))))

exp.df <- read_table("../data/meta/experimental-design.tsv", show_col_types = FALSE) %>% filter(sample %in% todo)
tmp <- exp.df %>%
    group_by(pool, group) %>%
    summarize(n.left = n())

write_csv(tmp, "todo-array.csv")