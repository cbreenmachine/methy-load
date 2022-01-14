library(data.table)
library(tidyverse)
library(wiscR)
library(parallel)
library(viridis)
options(scipen = 100)

ROOT_DIR <- "../../data/cov_meth/"
my_files <- list.files(path = ROOT_DIR, pattern =  "cov_meth.tsv", recursive = TRUE, full = TRUE)
# Make optional...?
my_files <- my_files[!str_detect(my_files, "X|Y")]

print(my_files)
print(paste0("Working on ", length(my_files), " total files"))

get_missing_counts <- function(file_name){
    # Loads in file; pivots it to be one position per row; one sample's coverage per column
    # Then returns a vector giving the tallies of missing values
    DT <- fread(file_name)
    DT.wide <- dcast(DT, sample~pos, value.var="coverage")
    pos.sorted <- sort(as.numeric(names(DT.wide)[-1]))

    # Set column order
    setcolorder(DT.wide, c("sample", as.character(pos.sorted)))
    pos <- as.numeric(names(DT.wide)[-1])
    N <- length(pos)
    
    # Difference: index 1 gives pos[2] - pos[1]
    dd <- pos[-1] - pos[1:(N-1)]

    left <- dd[1:(N-2)]
    right <- dd[2:(N-1)]
    min.dist <- pmin(left,right)

    num_nas <- colSums(is.na(DT.wide))[-1] # -1 gets rid of sample identifier (none of these will be na)
    # Trim off the first and last positions
    N <- length(num_nas)
    num_nas <- num_nas[2:(N-1)]

    out.df <- data.frame(count = as.vector(num_nas), min.dist = min.dist, left = left, right = right)
    return(out.df)
}

# Cobbles all missing values into one data frame (from each chromosome)
out <- mclapply(my_files, get_missing_counts, mc.cores = 6)
num_nas.df <- rbindlist(out)

# Simple labeling function to clean up long number strings
lab_fun <- function(x) paste0(x %/% 1000000, "MM")

# Some summary stats we want in the caption
N_CpGs <- lab_fun(nrow(num_nas.df))
N_CpGs_missing <- lab_fun(sum(num_nas.df$count > 0))

p <- num_nas.df %>%
    ggplot(aes(x = count)) +
    geom_histogram(fill = "#C5050C", binwidth = 4) +
    wiscR::light_theme() +
    scale_y_continuous(labels = lab_fun) +
    labs(caption = paste0("Number of CpGs: ", N_CpGs, 
    "\nNumber w/ at least one missing value: ", N_CpGs_missing)) +
    xlab("Number of samples with missing estimate") +
    ylab("Number of CpGs")


odir <- "./figs-missing/"
dir.create(odir, showWarnings = FALSE)

ofile <- file.path(odir, paste0(Sys.Date(), "-missing.png"))
wiscR::save_plot(p, ofile)

q <- num_nas.df %>%
    ggplot(aes(x = count, y = min.dist)) +
    geom_hex(aes(fill = stat((log(density))))) +
    scale_y_continuous(labels = function(x) paste0(x %/% 1000, "K")) +
    wiscR::light_theme() +
    scale_fill_viridis(option = "magma") +
    xlab("Number of samples with missing estimate") +
    ylab("Distance to closest neighbor (nt)")

ofile <- file.path(odir, paste0(Sys.Date(), "-dist-missing.png"))
wiscR::save_plot(q, ofile)