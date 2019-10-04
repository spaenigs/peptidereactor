.libPaths(c("apps/", .libPaths()))

library(yaml)
library(Interpol)

orig_list <- yaml.load_file(snakemake@input[["orig"]], as.named.list = FALSE)
enco_list <- yaml.load_file(snakemake@input[["enco"]], as.named.list = FALSE)

res <- Interpol(enco_list,
                dim = median(unlist(Map(nchar, orig_list))),
                method="spline")
df <- data.frame(res, row.names = attr(enco_list, "keys"))

write.csv(df, file = snakemake@output[[0]])