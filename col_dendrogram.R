suppressPackageStartupMessages(library(dplyr))
library(jsonlite)
library(data.tree)

df1 <- read.csv("heatmap_data.csv", row.names = 1, stringsAsFactors = FALSE)

hr <- hclust(dist(t(df1)), method = "average", members = NULL)
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

#df_res["group"] <- sapply(df_res[, "name"], function(x) {
#  if (length(grep("lambda-corr", x)) == 1 || length(grep("g-gap", x)) == 1) {
#    "psekraac"
#  } else if (nchar(x) == 0) {
#    "inner_node"
#  } else {
#    substr(x, 1, 6)
#  }
#})

lvls <- rownames(df1)
lvl_cnts <- seq_along(lvls)

colors <- mapply(function(x, y) c(x, y), lvls, lvl_cnts)

#for (l in lvls) {
#  df_res[df_res["group"] == l, "color"] <- colors[, l][2]
#}

df_res[df_res["name"] == "Root", "name"] <- ""

write_json(df_res, pretty = T, "dataset_col_out.json")