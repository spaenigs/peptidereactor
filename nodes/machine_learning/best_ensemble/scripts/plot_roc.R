library(yaml)
library(pROC)
  
y_true <- yaml.load_file(snakemake@input[["ytrue"]])
y_prob <- yaml.load_file(snakemake@input[["yprob"]])

png(snakemake@output[[1]])

rocobj <- plot.roc(y_true, y_prob,
                   main = "Confidence intervals", 
                   percent=TRUE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE)           # print the AUC (will contain the CI)
  
ciobj <- ci.se(rocobj,                         # CI of sensitivity
               specificities = seq(0, 100, 5)) # over a select set of specificities

plot(ciobj, type = "shape", col="lightgray", lty="dotted")  
plot(ci(rocobj, of = "thresholds", thresholds = "best"))

dev.off()

