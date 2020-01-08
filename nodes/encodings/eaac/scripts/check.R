library(yaml)

if (file.info(snakemake@input[["enco"]])$size == 0) {
  file.create(snakemake@output[[1]])
} else {
  source("nodes/encodings/eaac/scripts/interpolate.R")
}