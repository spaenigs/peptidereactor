library(MatrixCorrelation)

f <- Sys.glob("data/hiv_protease/csv/psekraac_type*/*.csv")
f <- f[1:20]
print(length(f))

m <- matrix(-1, nrow = length(f), ncol = length(f))
colnames(m) <- f
rownames(m) <- f

name_pairs <- expand.grid(x = f, y = f)

# TODO implement multi-threaded

for (i in 1:dim(name_pairs)[1]) {
  n1 <- as.character(name_pairs[i, "x"])
  n2 <- as.character(name_pairs[i, "y"])
  if (m[n1,n2] != -1) {
    #m[n2,n1] = m[n1,n2]
  } else if (m[n2,n1] != -1) {
    #m[n1,n2] = m[n2,n1]
  } else if (m[n1,n2] == -1) {
    X1 <- read.csv(n1, row.names = 1)
    X2 <- read.csv(n2, row.names = 1)
    m[n1,n2] <- 1 - RVadjMaye(data.matrix(X1[ , -which(names(X1) == "y")]),
                              data.matrix(X2[ ,-which(names(X2) == "y")]))
    m[n2,n1] <- m[n1, n2]
  }

}

hr <- hclust(dist(m), method = "complete", members=NULL)

#plot(hr)

par(mar=c(5,4,4,30))
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=1), horiz=T)