library(kaos)
library(yaml)

seq_yaml <- yaml.load_file(snakemake@input[[1]])
cgr_obj <- cgr(unlist(strsplit(seq_yaml$seq, "")), res=seq_yaml$res, sf=seq_yaml$sf)
seq_yaml$vec <- vectorize(cgr_obj)

write_yaml(seq_yaml, snakemake@output[[1]])

