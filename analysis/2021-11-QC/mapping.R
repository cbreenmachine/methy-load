library(tidyverse)
library(wiscR)

ROOT_DIR = "../../data/"
all_files <- list.files(ROOT_DIR, pattern = "[0-9][0-9][0-9].html", recursive = TRUE, full = TRUE)

get_mapping <- function(file){
    bn <- basename(file) %>% str_remove(".html")
    
    a <- readLines(file)
    ix <- grep("% Unique", a)
    mapping_perc <- str_remove(a[ix], "<TD>% Unique</TD><TD>") %>% str_remove("%</TD></TR>") %>% as.numeric()
    return(mapping_perc)
}

df <- data.frame(file = all_files) %>%
        dplyr::mutate(sample = str_remove( basename(file), ".html")) 
df$mapping <- as.numeric(sapply(all_files, get_mapping))

p <- df %>%
    ggplot(aes(x = mapping)) +
    geom_histogram(binwidth = 1) +
    wiscR::light_theme() +
    xlab("Unique mapping percentage") +
    ylab("Count")

wiscR::save_plot(p, "mapping.png")