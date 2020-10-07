library(foreach)
library(doParallel)
library(gtools)

cores <- detectCores()
cl <- makeCluster(cores)

registerDoParallel(cl)

c <- read.csv(snakemake@input[[1]], row.names = 1, stringsAsFactors = FALSE)

finalMatrix <- foreach(i=1:dim(c)[1], .combine=rbind) %dopar% {

  csv_path_1 <- c[i,1]
  csv_path_2 <- c[i,2]

  X1 <- read.csv(csv_path_1, row.names = 1)
  X2 <- read.csv(csv_path_2, row.names = 1)

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

  library(MatrixCorrelation)

  data.frame(
    e1 = csv_path_1,
    e2 = csv_path_2,
    R = RVadjMaye(m1, m2)
  )

}

stopCluster(cl)

write.csv(finalMatrix, snakemake@output[[1]])

