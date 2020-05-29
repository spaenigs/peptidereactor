library(scmamp)

d <- read.csv("res.csv", row.names = 1)

imanDavenportTest(d)

nm <- nemenyiTest(d)
nm$diff.matrix

svg(filename="cd.svg")
plotCD(results.matrix = d, alpha = 0.05)
dev.off()

test.res <- postHocTest(data = d, test ='friedman', correct ="holland")
test.res