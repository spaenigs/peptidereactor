suppressPackageStartupMessages(library(dplyr))
library(data.tree)
library(jsonlite)

df1 <- read.csv("rv.csv", row.names = 1, stringsAsFactors = FALSE) #[1:100, ]

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

df_res["name"] <- Map(trimws, df_res["name"])
df_res["name"] <- Map(function(x) gsub("[°-]", "", x), df_res["name"])
df_res["name"] <- Map(function(x) gsub(" ", "", x), df_res["name"])
df_res["name"] <- Map(function(x) gsub("¦", "", x), df_res["name"])
df_res["name"] <- Map(function(x) gsub("^(\\d+)", "", x), df_res["name"])

df_res["group"] <- sapply(as.vector(df_res[, "name"]), function(x) {
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

write_json(df_res, pretty = T, "rv.json")

