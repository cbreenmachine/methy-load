
# System data, make sankey diagram
# Make a python crawler that creates the extracted file if it's not there already
max_possible <- 100 

root_dir <- "../../data/"
raw_fastqs <- list.files(root_dir, pattern = "fastq.gz", recursive = TRUE)
trimmed_fastqs <- list.files(root_dir, pattern = ".fq", recursive = TRUE)
mapped <- list.files(root_dir, pattern = ".bam.csi", recursive = TRUE)
called <- list.files(root_dir, pattern = "^[0-9][0-9][0-9].bcf.csi", recursive = TRUE)
extracted <- list.files(root_dir, pattern = ".tsv", recursive = TRUE) 

N <- min(c(max_possible, length(raw_fastqs)))

df <- data.frame(NumTrimmed = length(trimmed_fastqs) / 2,
                 NumMapped  = length(mapped),
                 NumCalled = length(called),
                 NumExtracted = length(extracted))
print(df)
