suppressPackageStartupMessages(library(dplyr))
library(gtools)
library(yaml)
library(MatrixCorrelation)

res <- read_yaml(snakemake@input[[1]])

dfs <- lapply(res, function(e) {
  data.frame(x = e[1], y = e[2])
})

name_pairs <- bind_rows(dfs)

df_res <- data.frame()

for (i in 1:dim(name_pairs)[1]) {

    n1 <- as.character(name_pairs[i, 1])
    n2 <- as.character(name_pairs[i, 2])

    X1 <- read.csv(n1, row.names = 1)
    X2 <- read.csv(n2, row.names = 1)

    if (dim(X1)[1] > dim(X2)[1]) {
      m1 <- data.matrix(X1[rownames(X2), -which(names(X1) == "y")])
      m2 <- data.matrix(X2[            , -which(names(X2) == "y")])
    } else if (dim(X1)[1] < dim(X2)[1]) {
      m1 <- data.matrix(X1[            , -which(names(X1) == "y")])
      m2 <- data.matrix(X2[rownames(X1), -which(names(X2) == "y")])
    } else {
      m1 <- data.matrix(X1[ , -which(names(X1) == "y")])
      m2 <- data.matrix(X2[ , -which(names(X2) == "y")])
    }

    df <- data.frame(
      x = gsub(".*/(.*)\\.csv", "\\1", n1),
      y = gsub(".*/(.*)\\.csv", "\\1", n2),
      res = 1 - RVadjMaye(m1, m2)
    )

    df_res <- rbind(df_res, na.omit(df))
}

write.csv(df_res, snakemake@output[[1]])