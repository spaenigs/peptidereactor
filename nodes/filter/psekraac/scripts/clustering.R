library(cluster)

d <- read.csv(snakemake@input[[1]], row.names = "e1")

# get representative file by computing mediod of one cluster
km <- pam(d, 1)

write.csv(data.frame(row.names(km$medoids)), file = snakemake@output[[2]])
    
