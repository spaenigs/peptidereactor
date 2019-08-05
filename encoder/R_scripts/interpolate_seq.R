#!/usr/bin/env Rscript

#suppressWarnings(library(Interpol))
library(Interpol)

input <- commandArgs(trailingOnly = TRUE)
filePath = input[1]
interpolateTo <- as.numeric(input[2])

conn <- file(filePath,open="r")
inputSequences <- eval(parse(text=readLines(conn)[1]))
close(conn)

res = Interpol(inputSequences, dims = interpolateTo, method = "spline")
res_as_str <- paste0("[", paste(apply(res, 1, function(row) paste0("[", paste(row, collapse=", "), "]")), collapse=", "), "]")

cat(res_as_str)


