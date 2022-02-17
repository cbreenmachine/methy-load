library(tidyverse)

df <- read_tsv("../data/meta/ADRC-samplesheet-cleaned.tsv")

df

list.files("../data/", recursive = TRUE, full=TRUE, pattern = "04-extract/")

extracted <- list.files("../data/", pattern = "^[0-9][0-9][0-9].tsv", recursive = TRUE) 

todo <- as.character(sort(setdiff(df$ALISCH_ID, str_remove(basename(extracted), ".tsv"))))

exp.df <- read_table("../data/meta/experimental-design.tsv", sep = " ") 

exp.df


