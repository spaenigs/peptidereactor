library(yaml)
library(pROC)

#setClass("Foo", representation(input="list", output="list"))
#
#snakemake = new("Foo",
#                input=list(ytrue=c(0.0, 0.0), yprob=c(0.4744404761904762, 0.2406666666666667)),
#                output=list("roc.png"))
#
#y_true <- snakemake@input[["ytrue"]]
#y_prob <- snakemake@input[["yprob"]]

df1 = read.csv("data/hiv_protease/benchmark/ensemble/y_true_cv_egaac_window_8.csv", row.names = 1)
df2 = read.csv("data/hiv_protease/benchmark/ensemble/y_prob_cv_egaac_window_8.csv", row.names = 1)

y_true = as.vector(t(df1))
y_prob = as.vector(t(df2))

#png(snakemake@output[[1]])

svg("roc.svg")

rocobj <- plot.roc(y_true, y_prob,
                   main = "Confidence intervals",
                   percent=TRUE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE)           # print the AUC (will contain the CI)

ciobj <- ci.se(rocobj,                         # CI of sensitivity
               specificities = seq(0, 100, 5)) # over a select set of specificities

plot(ciobj, type = "shape", col="lightgray", lty="dotted")
plot(ci(rocobj, of = "thresholds", thresholds = "best"))

#plot(c(0, 5, 10, 15 , 20, 25 , 30 , 35 , 40 , 45,  50, 55 , 60 , 65, 70  ,75 , 80 , 85,  90,  95 , 100), rev(ciobj[ ,1]), type = "l")
#lines(c(0, 5, 10, 15 , 20, 25 , 30 , 35 , 40 , 45,  50, 55 , 60 , 65, 70  ,75 , 80 , 85,  90,  95 , 100), rev(ciobj[ ,2]))
#lines(c(0, 5, 10, 15 , 20, 25 , 30 , 35 , 40 , 45,  50, 55 , 60 , 65, 70  ,75 , 80 , 85,  90,  95 , 100), rev(ciobj[ ,3]))

dev.off()

