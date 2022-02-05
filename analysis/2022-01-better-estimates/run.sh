#Rscript estimate.R --sample --covariates mean_methylation,CD8T,CD4T,NK,Bcell,Mono,Gran --ofile ./data/PC-cell-compositions.csv
#Rscript estimate.R --sample --covariates mean_methylation,CD8T,CD4T,NK,Bcell,Mono,Gran,age,sex,bmi --ofile ./data/PC-cell-compositions-pheno.csv

Rscript plot-estimate-PCs.R --ifile ./data/PC-cell-compositions.csv --odir ./figs/PC-cell-compositions/
Rscript plot-estimate-PCs.R --ifile ./data/PC-cell-compositions-pheno.csv --odir ./figs/PC-cell-compositions-pheno/