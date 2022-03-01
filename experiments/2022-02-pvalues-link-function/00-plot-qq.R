library(tidyverse)
library(wiscR)
library(latex2exp)

ifile <- "../../analysis/2022-01-diff-analysis/chr18/smooth-100-PCs-2/models.RData"
load(ifile)

y <- test.cohort$pvals
y <- y[!is.na(y)]

p.obs <- quantile(y, probs = seq(0,1, length = length(y)))
x <- seq(0,1, length=length(y))

p <- data.frame(x = -log10(x), y = -log10(p.obs)) %>%
    ggplot(aes(x = x, y= y)) +
    geom_point(size=2.5) +
    wiscR::light_theme() +
    geom_abline(intercept = 0, slope = 1, color = "blue", size = 2) +
    xlab(TeX("Expected quantile of $-log_{10}(p)$ if $P \\sim Unif(0,1)$")) +
    ylab(TeX("Observed quantile of $-log_{10}(p)$"))
wiscR::save_plot(p, "chr18-smooth100-pcs2-qqunif.png")


p <- data.frame(x = (x), y = (y)) %>%
    ggplot(aes(x = x, y= y)) +
    geom_point(size=2.5) +
    wiscR::light_theme() +
    geom_abline(intercept = 0, slope = 1, color = "blue", size = 2) +
    xlab(TeX("Expected quantile of $-log_{10}(p)$ if $P \\sim Unif(0,1)$")) +
    ylab(TeX("Observed quantile of $-log_{10}(p)$"))
wiscR::save_plot(p, "tmp.png")


p <- data.frame(y = y) %>%
    ggplot(aes(x = y, y = ..density..)) +
    geom_histogram(bins = 50) +
    wiscR::light_theme() +
    xlab("Raw P-value") +
    ylab("Density") +
    geom_hline(yintercept = 1, color = "red", size = 2)
wiscR::save_plot(p, "chr18-smooth100-pcs2-hist.png")


png("qqnorm.png")
stats::qqnorm(-log10(y))
dev.off()