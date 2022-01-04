library(data.table)
library(tidyverse)
library(wiscR)
library(parallel)
options(scipen = 100)

ROOT_DIR <- "../../data/2021-12-02-methylation-coverage/"
my_files <- list.files(path = ROOT_DIR, pattern =  "meth_cov.tsv", recursive = TRUE, full = TRUE)

get_missing_counts <- function(file_name){
    DT <- fread(file_name)
    DT.wide <- dcast(DT, sample~pos, value.var="coverage")

    num_nas <- colSums(is.na(DT.wide))[-1] # -1 gets rid of sample identifier (none of these will be na)
    return(as.vector(num_nas))
}

out <- mclapply(my_files, get_missing_counts, mc.cores = 6)
num_nas.df <- data.frame(count = unlist(out))

lab_fun <- function(x) round(x, -5)

N_CpGs <- round(nrow(num_nas.df), -5)
N_CpGs_missing <- round(sum(num_nas.df$count > 0), -5)

p <- num_nas.df %>%
    ggplot(aes(x = count)) +
    geom_bar(fill = "#C5050C") +
    wiscR::light_theme() +
    scale_y_continuous(labels = lab_fun) +
    labs(caption = paste0("Number of CpGs: ", N_CpGs, 
    "\nNumber w/ at least one missing value: ", N_CpGs_missing)) +
    xlab("Number of samples with missing estimate") +
    ylab("Number of CpGs")


wiscR::save_plot(p, "missing.png")