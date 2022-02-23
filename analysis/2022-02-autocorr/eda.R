
library(data.table) 
library(tseries)
library(Rlab)
library(wiscR)
library(tidyverse)
# Questions: what is the autocorrelation here????
# What about with p-values??


DT <- fread("../2022-01-diff-analysis/chr18/methylation.csv")
xx <- as.vector(unlist(DT[ , "180", with=F]))
xx <- xx[!is.na(xx)]
ac.xx <- as.vector(acf(xx, pl=FALSE, lag = 1000)[[1]])

#--> What baout null
yy <- rbern(n = length(xx), prob = mean(xx))
ac.yy <- as.vector(acf(yy, pl=FALSE)[[1]])


get_ac <- function(ix){
    xx <- as.vector(unlist(DT[ , ix, with=F]))
    xx <- xx[!is.na(xx)]
    ac.xx <- as.vector(acf(xx, pl=FALSE, lag = 1000)[[1]])
    return(ac.xx)
}

data <- lapply(2:10, get_ac)
df <- as.data.frame(do.call("cbind", data))
names(df) <- names(DT)[2:10]
df$lag <- 0:(nrow(df)-1)

#--> Plot real data 
p <- df %>%
    pivot_longer(cols = -lag, values_to = "auto.corr", names_to = "sample") %>%
    ggplot(aes(x = lag, y = auto.corr, color =sample)) +
    geom_smooth(se=F) +
    wiscR::light_theme()

wiscR::save_plot(p, "autocorr.png")

