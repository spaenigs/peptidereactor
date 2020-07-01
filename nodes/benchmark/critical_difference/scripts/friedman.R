library(scmamp)
library(yaml)

Sys.setlocale("LC_NUMERIC","en_US.UTF-8")

d <- read.csv(snakemake@input[[1]], row.names = 1)

idt <- imanDavenportTest(d)

write_yaml(
  list(
    method = idt$method, statistic = idt$statistic, p.value = idt$p.value
  ), snakemake@output[["fo"]]
)

nm <- nemenyiTest(d)
write_yaml(
  list(
    method = nm$method, cd = nm$statistic
  ), snakemake@output[["no"]]
)

dm <- nm$diff.matrix
rownames(dm) <- colnames(dm)
write.csv(dm, snakemake@output[["co"]])

#svg(filename = "cd.svg")
# algorithms that cannot be regarded as different are joined by a horizontal bar
#plotCD(results.matrix = d, alpha = 0.05)
#dev.off()

#test.res <- postHocTest(data = d, test = 'friedman', correct = "holland")
#test.res

