library(data.table)
library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-odir", default=paste0("figs-", Sys.Date(), "/"))
parser$add_argument("-chrom", default="6", help="Which chromosome to run analysis on.")
parser$add_argument("-mc.dir", default="../../data/cov_meth/")
parser$add_argument("-start", default=3.3277e7)
parser$add_argument("-stop", default=3.3284e7)
args <- parser$parse_args()
dir.create(args$odir, showWarning = FALSE)

ifile <- file.path(args$mc.dir, paste0("chr", args$chrom, "_cov_meth.tsv"))

DT <- fread(ifile)
DT.sub <- DT[args$start < pos]
DT.sub <- DT.sub[pos < args$stop]
DT.sub$methylation <- DT.sub$methylated / DT.sub$coverage

DT.wide <- dcast(DT.sub, sample~pos, value.var="methylation")
print(nrow(DT.sub))


df <- read_csv("../../data/meta/phenos-cleaned.csv")

p <- DT.sub %>%
    inner_join(df, by = "sample") %>%
    ggplot(aes(x = pos, y = methylation, group = sample, color = cohort)) +
    geom_smooth(se = F) +
    geom_point(alpha = 0.3) 

ggsave("tmp.png", p)
