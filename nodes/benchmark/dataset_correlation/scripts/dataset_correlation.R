library(gtools)
library(foreach)
library(doParallel)

files <- c(
  Sys.glob(paste0(snakemake@input[["group_1"]], "*.csv")),
  Sys.glob(paste0(snakemake@input[["group_2"]], "*.csv"))
)

name_pairs <- combinations(length(files),2, files, repeats=FALSE)

cl <- makeCluster(snakemake@params$cores)
registerDoParallel(cl)

chunks <- split(1:dim(name_pairs)[1], sort(1:dim(name_pairs)[1]%%10000))

df_res <- data.frame()

for (chunk in chunks) {

  finalMatrix <- foreach(i=chunk, .combine=rbind) %dopar% {

    library(MatrixCorrelation)

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

    na.omit(df)
  }

  token <- paste0(sample(c(0:9, LETTERS[1:6]), 8, T), collapse = '')
  df_res <- rbind(df_res, finalMatrix)
  write.csv(df_res, paste0(snakemake@output[[2]], "/dataset_correlation_part_", token, ".csv"))
}

stopCluster(cl)

write.csv(df_res, snakemake@output[[1]])