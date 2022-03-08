suppressPackageStartupMessages({
    library(viridis)
    library(tidyverse)
    library(data.table)
    library(wiscR)
    library(DSS)
    library(fdrtool)

}) 


# Load data --------------------------------------------------------------------
all.files <- system("find . -name models.RData | grep smooth-150-PCs-2", intern = T)

.correct_pvals <- function(df){
    zz <- (df$stat - mean(df$stat)) / sd(df$stat)
    out <- fdrtool(zz, statistic = "normal", plot = FALSE, cutoff.method = "locfdr")

    df$pvals.adj <- out$pval
    df$lfdr <- out$lfdr 
    df$qval <- out$qval
    return(df)
}


.load_data <- function(path){
    # Load from path like "./chr8/smooth-150-PCs-2/models.RData"
    load(path)
    df <- data.frame(chr = test.cohort$chr, pos = test.cohort$pos, 
                stat = test.cohort$stat) %>% drop_na() 
    return(df)
}

.wrapper <- function(path){ .load_data(path) %>% .correct_pvals() %>% return()}
df <- do.call("rbind", lapply(all.files, .wrapper))  %>% mutate(p = pvals.adj)



# Plot functions -------------------------------------------------------------------
# Thanks to https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
.pad_chrom <- function(s="chr1"){
    # chr1 --> chr01
    t <- str_remove(s, "chr") %>% 
        str_pad(width = 2, pad = "0") 
    paste0("chr", t) %>% return()
}


.thin_data <- function(df, alpha=0.1, prop=0.05){
    # Samples a proportion of non-significant sites to make plotting faster
    sig.df <- df %>% dplyr::filter(p < alpha)
    nonsig.df <- df %>% 
        dplyr::filter(p > alpha) %>% 
        dplyr::slice_sample(prop = prop)
    
    thinned.df <- rbind(sig.df, nonsig.df) %>% 
        mutate(chr = .pad_chrom(chr)) 
    return(thinned.df)
}


.adjust_pos <- function(df){
    data.cum <- df %>% 
        group_by(chr) %>% 
        summarise(max_bp = max(pos)) %>% 
        mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
        dplyr::select(chr, bp_add)

    
    df %>% 
        inner_join(data.cum, by = "chr") %>% 
        mutate(bp_cum = pos + bp_add) %>% return()
}

.generate_axis_set <- function(df){
    df %>% group_by(chr) %>% 
        summarize(center = mean(bp_cum)) %>% return()
}

sub.df <- .thin_data(df) %>% .adjust_pos()
axis.set <- .generate_axis_set(sub.df)

p <- sub.df %>% 
  ggplot(aes(x = bp_cum, color = as.factor(chr), y = -log10(p))) +
  geom_point(alpha = 0.9, size = 3.5) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  wiscR::light_theme() +
  xlab("") +
  ylab("-log10(p)") +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, size = 20, vjust = 0.5)) +  
  geom_hline(yintercept = 0,  size = 1, color = "grey") +
  theme(plot.caption = element_text(hjust = 0, face= "italic", size = 20),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA)) +
  scale_color_manual(values = rep(c("grey", "#C5050C"), 12)) +
  xlab("Genomic position")

wiscR::save_plot(p, "manhattan-adjusted.png")
  