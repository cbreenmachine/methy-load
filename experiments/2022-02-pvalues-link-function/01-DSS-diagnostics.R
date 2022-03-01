# 01-find-DMRs.R
# Runs Dispersion Shrinkage Estimation method on methylation data (one chromosome at a time)
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(argparse)
    library(DSS)
    library(fastglm)
    library(wiscR)
})

parser <- ArgumentParser()
# Things that may will change
parser$add_argument("--chr", default= "chr18", help='Chromosome to run DSS on')
parser$add_argument("--covariates", default= "Gran,CD8T,CD4T,NK,Bcell,bmi,age,sex,PC1,PC2")
parser$add_argument("--smoothing", default= 100, help= 'Width of smoothing window')
# File paths
parser$add_argument("--samples_file", default= "../../data/meta/array-samples.csv", help="CSV file with samples to filter")
parser$add_argument("--pheno_file", default= "../../data/meta/phenos-cleaned.csv", help="CSV file with LOAD/control status")
parser$add_argument("--pc_dir", default= "../../data/prin-comps-array-samples/", help="CSV file with LOAD/control status")
args <- parser$parse_args()

n.sample <- 1000


#--> Directories
idir <- file.path("../../analysis/2022-01-diff-analysis/", args$chr)
odir <- "./"
dir.create(odir, showWarnings = FALSE)

#--> Handle Models
covariates <- as.vector(unlist(str_split(args$covariates, ",")))
mod.form <- as.formula(paste(c("~cohort", covariates), collapse = "+"))

#--> Read data
pc.df <- read_csv(file.path(args$pc_dir,  paste0(args$chr, ".csv")), col_types = cols())
ss.df <- read_csv(args$pheno_file, col_types = cols())
M <- fread(file.path(idir, "M.csv"))
Cov <- fread(file.path(idir, "Cov.csv"))

#--> Check that we have phenotypes/PCs for all columns in M/Cov
valid.samples <- intersect(intersect(colnames(M), ss.df$sample), pc.df$sample)
filt.df <- ss.df %>%
            inner_join(pc.df, by = "sample") %>% 
            dplyr::filter(sample %in% valid.samples) 
drop.cols <- as.character(setdiff(colnames(M), valid.samples))[-1] # dont remove pos

if (length(drop.cols) > 0){
    M[, (drop.cols):=NULL]
    Cov[, (drop.cols):=NULL]
}

keep.ix <- which(rowSums(Cov == 0) ==  0)
M <- M[keep.ix, ]
Cov <- Cov[keep.ix, ]

#--> Check ordering
all(filt.df$sample == names(M)[-1])

#--> Sampling
set.seed(919)
sample.ix <- sample(1:nrow(M), size = n.sample, replace = FALSE)

M.samp <- M[sample.ix, ]
Cov.samp <- Cov[sample.ix, ]

# create bs seq object, needs chromosome identiifer, methylated reads, and unmethylated reads
bs <- BSseq(chr = rep(args$chr, nrow(M.samp)), pos = M.samp$pos,
            M = as.matrix(M.samp[ , -c("pos"), with=FALSE]), 
            Cov = as.matrix(Cov.samp[, -c("pos"), with=FALSE]), 
            sampleNames = names(M.samp)[-1])
all( filt.df$sample == colnames(bs) )


filt.df$cohort <- as.factor(filt.df$cohort)
dml.fit <- DMLfit.multiFactor(bs, design = filt.df, smoothing = args$smooth, 
            smoothing.span = args$smoothing, formula = mod.form)

colnames(dml.fit$X)[2]
beta.LOAD.test <- dml.fit$fit$beta[ ,2]


#--> Plot residuals
observed <- M.samp[ , -"pos", with=FALSE] / Cov.samp[ , -"pos", with=FALSE]

get_res <- function(ix){
    predicted <- dml.fit$X %*% dml.fit$fit$beta[ix, ]
    res <- unlist(observed[ix, ] - predicted)
    # sq.res <- res ^2
    return(res)
}

res.df <- as.data.frame(do.call(cbind, lapply(1:1000, get_res))) %>% mutate(method = "DSS")
dir.create("residuals/")
plot_by_ix <- function(ix){
    p <- res.df %>%
        ggplot(aes_string(x=names(.)[ix])) +
        geom_histogram() +
        wiscR::light_theme() +
        xlim(c(-1, 1)) +
        facet_wrap(.~ method)
    wiscR::save_plot(p, file.path("residuals/", paste0("V", ix, ".png")))
}

lapply(1:20, plot_by_ix)

# permute_DSS <- function(i){
#     tmp.fit <- DMLfit.multiFactor(bs, design = mutate(filt.df, cohort = sample(cohort)), smoothing = args$smooth, 
#                 smoothing.span = args$smoothing, formula = mod.form)
#     print(colnames(tmp.fit$X)[2])
#     beta <- tmp.fit$fit$beta[ ,2]
#     print(i)
#     return(beta)
# }

# permutations <- lapply(1:1000, permute_DSS)
# perm.df <- as.data.frame(do.call(rbind, permutations))


# write_csv(perm.df, "./permutations.beta.csv")
# write_csv((as.data.frame(beta.LOAD.test)), "./test.beta.csv")


# get_p <- function(ix){
#     N <- sum(abs(perm.df[ ,ix])  > abs(beta.LOAD.test[ix])) + 1
#     D <- nrow(perm.df) + 1
#     return( N / D)
# }

# p.vals <- as.data.frame(do.call(rbind, lapply(1:1000, get_p)))
# write_csv(p.vals, "p.vals-permutation-DSS.csv")


# p <- p.vals %>%
#     ggplot(aes(x = V1, y = ..density..)) +
#     geom_histogram() +
#     wiscR::light_theme() +
#     xlab("P-value") +
#     ylab("Density") +
#     labs(caption = "Permutation p-values on 1000 CpGs on chr 18")

# wiscR::save_plot(p, "p.vals-permutation-DSS.png")


fit_logit <- function(ix){
    x <- model.matrix( mod.form, data = filt.df)
    y <- unlist(M.samp[ix, -"pos", with=F] / Cov.samp[ix, -"pos", with = F])
    
    glm.fitted <- fastglm::fastglm(y = y, x = x, family = quasibinomial(link = "logit"))
    p <- coefficients(summary(glm.fitted))[2 ,4]
    res <- y - predict(glm.fitted, x)
    return(list(p, res))
}

out <- lapply(1:1000, fit_logit)

logit.p <-  as.data.frame(do.call(rbind, lapply(out, `[[`, 1)))
res.logit.df <-  as.data.frame(do.call(cbind, lapply(out, `[[`, 2))) %>% mutate(method = "logit")

res.both.df <- cbind(res.df, res.logit.df)


write_csv(logit.p, "p.vals-logit-link.csv")

p <- logit.p %>%
    ggplot(aes(x = V1, y = ..density..)) +
    geom_histogram() +
    wiscR::light_theme() +
    xlab("P-value") +
    ylab("Density") +
    labs(caption = "Logit link p-values on 1000 CpGs on chr 18")

wiscR::save_plot(p,  "p.vals-logit-link.png")
