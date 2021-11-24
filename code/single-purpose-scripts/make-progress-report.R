
# System data, make sankey diagram
# Make a python crawler that creates the extracted file if it's not there already
max_possible <- 100 

root_dir <- "../../data/"
raw_fastqs <- list.files(root_dir, pattern = "fastq.gz", recursive = TRUE)
trimmed_fastqs <- list.files(root_dir, pattern = ".fq", recursive = TRUE)
mapped <- list.files(root_dir, pattern = ".bam.csi", recursive = TRUE)
extracted <- list.files(root_dir, pattern = ".tsv", recursive = TRUE) 

N <- min(c(max_possible, length(raw_fastqs)))

df <- data.frame(PercTrimmed = length(trimmed_fastqs) / (2*N),
                 PercMapped  = length(mapped) / N,
                 PercExtracted = length(extracted) / N)
print(df)
