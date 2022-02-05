# Want a matrix that is chromosome
# https://stats.stackexchange.com/questions/24934/command-line-tool-to-calculate-basic-statistics-for-stream-of-values
# Use alias
bcftools query --format '[ %DP]' --include 'CG="Y"' 101_chr22.bcf | summary