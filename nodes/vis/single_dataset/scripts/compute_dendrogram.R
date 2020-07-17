suppressPackageStartupMessages(library(dplyr))
library(data.tree)
library(jsonlite)

df_f1 <- read.csv(snakemake@input[[1]], row.names = 1, check.names = F)

df_groups <- data.frame(encoding = colnames(df_f1))
df_groups["group"] <- sapply(as.vector(df_groups[, "encoding"]), function(x) {
  if (length(grep("lambda-corr", x)) == 1 || length(grep("g-gap", x)) == 1) {
    "psekraac"
  } else if (nchar(x) == 0) {
    "inner_node"
  } else {
    substr(x, 1, 6)
  }
})

df_groups["median"] <- apply(df_f1, 2, median)
df_groups <- df_groups[order(-df_groups["median"]), ]
top_encodings <- df_groups[1:50, "encoding"]

df1 <- read.csv(snakemake@input[[2]], row.names = 1, stringsAsFactors = FALSE)
df1 <- df1 %>% filter(x %in% top_encodings) %>% filter(y %in% top_encodings)

encodings <- unique(c(df1$x, df1$y))

m <- matrix(-1, nrow = length(encodings), ncol = length(encodings))
colnames(m) <- encodings
rownames(m) <- encodings

for (i in 1:dim(df1)[1]) {
  x <- df1[i, "x"]
  y <- df1[i, "y"]
  m[x, y] <- df1[i, "res"]
  m[y, x] <- df1[i, "res"]
}

hr <- hclust(dist(m), method = "complete", members = NULL)
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

df_res["group"] <- sapply(df_res[, "name"], function(x) {
  if (length(grep("lambda-corr", x)) == 1 || length(grep("g-gap", x)) == 1) {
    "psekraac"
  } else if (nchar(x) == 0) {
    "inner_node"
  } else {
    substr(x, 1, 6)
  }
})

lvls <- unique((df_res %>%
  filter(group != "inner_node") %>%
  filter(group != "Root"))[, "group"])
lvl_cnts <- seq_along(lvls)

colors <- mapply(function(x, y) c(x, y), lvls, lvl_cnts)

for (l in lvls) {
  df_res[df_res["group"] == l, "color"] <- colors[, l][2]
}

df_res[df_res["name"] == "Root", "name"] <- ""

write_json(df_res, pretty = T, snakemake@output[[1]])

