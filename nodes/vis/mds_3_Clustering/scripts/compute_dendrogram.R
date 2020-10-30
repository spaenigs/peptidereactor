suppressPackageStartupMessages(library(dplyr))
library(jsonlite)
library(data.tree)

df1 <- read.csv(snakemake@input[[1]], row.names = 1, stringsAsFactors = FALSE)

if (snakemake@wildcards[["axis"]] == "row") {
  dist_data <- dist(df1)
} else {
  dist_data <- dist(t(df1))
}

hr <- hclust(dist_data, method = "average", members = NULL)
de <- as.dendrogram(hr)
n <- as.Node(de)

n$Do(function(node) node$nameId <- paste0(sample(c(0:9, LETTERS[1:6]), 8, T), collapse = ''))
n$Do(function(node) node$parentId <- node$parent$nameId)

df_res <- ToDataFrameTree(n, "nameId", "parentId")
df_res["id"] <- rownames(df_res)

for (i in rownames(df_res)) {
  if (!is.na(df_res[i, "parentId"])) {
    childOf <- df_res[trimws(df_res[i, "parentId"]) == trimws(df_res[, "nameId"]), "id"]
    df_res[i, "childOf"] <- childOf
  }
}

colnames(df_res) <- c("name", "nameId", "parentId", "id", "parent")
df_res <- df_res[, c("name", "id", "parent")]

df_res["name"] <- sapply(df_res[ ,"name"], function(x) {
  x <- trimws(x)
  x <- gsub("(.*?)(\\w+.*\\w+)(.*?)", "\\2", x, perl = T)
  x <- gsub("[¦°]", "", x)
  x <- gsub("--\\d*", "", x)
  gsub(" ", "", x)
})

lvls <- rownames(df1)
lvl_cnts <- seq_along(lvls)

colors <- mapply(function(x, y) c(x, y), lvls, lvl_cnts)

df_res[df_res["name"] == "Root", "name"] <- ""

if (snakemake@wildcards[["axis"]] == "col") {
  name <- "tree_col"
} else {
  name <- "tree"
}

write_json(df_res, pretty = T, snakemake@output[[1]])

#json <-  toJSON(list(name=name, values=df_res), pretty = TRUE)
#gs <- gsub(paste0("[\"", name, "\"]"), paste0("\"", name, "\""), json, fixed=TRUE)
#
#file <- file(snakemake@output[[1]])
#writeLines(gs, file)
#close(file)