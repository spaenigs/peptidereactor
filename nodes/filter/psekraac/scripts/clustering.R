library(cluster)

d <- read.csv(snakemake@input[[1]], row.names = "e1")

# get representative file by computing mediod if one cluster
km <- pam(d, 1)

png(snakemake@output[[1]], width=800, height = 500)

plot(d[ , "x1"], d[ , "x2"], pch=4, main = basename(row.names(km$medoids)), xlab = "tSNE 1", ylab = "tSNE 2")
points(x=km$medoids[1, 1], y=km$medoids[1, 2], col="blue", pch=19)

dev.off()

write.csv(data.frame(row.names(km$medoids)), file = snakemake@output[[2]])
    
