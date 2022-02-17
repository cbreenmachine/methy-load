library(tidyverse)
library(tools)

avail.files <- system("ls ../data/batch01/pool0[1-5]-group0[1-5]/04-extract/[0-9][0-9][0-9].tsv", intern = TRUE)
avail.samples <- file_path_sans_ext(basename(avail.files))

df <- read_csv("../data/meta/phenos-cleaned.csv", show_col_types = FALSE) %>%
        filter(sample %in% avail.samples)

