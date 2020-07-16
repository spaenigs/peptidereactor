library(gtools)
library(foreach)
library(doParallel)

f <- c(
  Sys.glob("data/hiv_protease/csv/sequence_based/*.csv"),
  Sys.glob("data/hiv_protease/csv/structure_based/*.csv")
)
#f <- f[1:5]

#print(length(f))

#m <- matrix(-1, nrow = length(f), ncol = length(f))
#colnames(m) <- f
#rownames(m) <- f

#name_pairs <- expand.grid(x = f, y = f)
name_pairs <- combinations(length(f),2, f, repeats=FALSE)

cl <- makeCluster(10)
registerDoParallel(cl)

for (chunk in split(1:dim(name_pairs)[1], sort(1:dim(name_pairs)[1]%%10000))) {
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
   df_tmp <- read.csv("rv.csv", row.names = 1)
   write.csv(rbind(df_tmp, finalMatrix), "rv.csv")
}




stopCluster(cl)


# df1 = finalMatrix
#m <- matrix(-1, nrow = length(unique(c(df1$x, df1$y))), ncol = length(unique(c(df1$x, df1$y))))
#colnames(m) <- unique(c(df1$x, df1$y))
#rownames(m) <- unique(c(df1$x, df1$y))

#for (i in 1:dim(df1)[1]) {
#  x = df1[i, "x"]
#  y = df1[i, "y"]
#  m[x,y] = df1[i, "res"]
#  m[y,x] = df1[i, "res"]
#}

#for (i in 1:dim(name_pairs)[1]) {
#  n1 <- as.character(name_pairs[i, "x"])
#  n2 <- as.character(name_pairs[i, "y"])
#  if (m[n1,n2] != -1) {
#    #m[n2,n1] = m[n1,n2]
#  } else if (m[n2,n1] != -1) {
#    #m[n1,n2] = m[n2,n1]
#  } else if (m[n1,n2] == -1) {
#    X1 <- read.csv(n1, row.names = 1)
#    X2 <- read.csv(n2, row.names = 1)
#    m[n1,n2] <- 1 - RVadjMaye(data.matrix(X1[ , -which(names(X1) == "y")]),
#                              data.matrix(X2[ ,-which(names(X2) == "y")]))
#    m[n2,n1] <- m[n1, n2]
#  }
#
#}
#
#hr <- hclust(dist(m), method = "complete", members=NULL)
# de = as.dendrogram(hr)
##plot(hr)
#
#par(mar=c(5,4,4,30))
#plot(as.dendrogram(hr), edgePar=list(col=3, lwd=1), horiz=T)
#
# attributes(de[[2]])
# attr(de[[1]], "members")