library(tidyverse)

df <- read_csv("../../data/meta/phenos-cleaned-all.csv")

df.2 <- df %>% 
        select(c(CD8T, CD4T, NK, Bcell, Mono, Gran, ab_142, p_tau, h_tau,
                age, pack_years, charlson_score, rey15_raw))

df.2 <- df %>% 
        select(c(CD8T, CD4T, NK, Bcell, Mono, Gran, 
                age, charlson_score, rey15_raw, bmi)) %>%
                as.data.frame()
rownames(df.2) <- df$Alisch


ix <- is.na(df.2$rey15_raw)
df.2$rey15_raw[ix] <- floor(mean(df.2$rey15_raw[!ix]))

for (vv in colnames(df.2)){
    xx <- df.2[[vv]]
    z <- (xx - mean(xx)) / sd(xx)
    print(vv)
    print( rownames(df.2)[which(abs(z) > 1.5)] )
}
