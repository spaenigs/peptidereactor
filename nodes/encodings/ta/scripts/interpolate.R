library(yaml)
library(Interpol)

enco_list <- yaml.load_file(snakemake@input[["enco"]])

res <- Interpol(enco_list[["enco_seqs"]],
                dim = enco_list[["interpolate_to"]],
                method="spline")
df <- data.frame(res, row.names = names(enco_list$enco_seqs))

write.csv(df, file = snakemake@output[[1]])