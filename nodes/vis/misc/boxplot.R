library(reshape2)
library(ggplot2)

df = read.csv("data/hiv_protease/benchmark/metrics/f1.csv", row.names = 1)
df["fold"] = 1:dim(df)[1]

df_melted = melt(df,
                 id.vars = "fold",
                 measure.vars = colnames(df)[1:dim(df)[2]-1],
                 variable.name = "Encoding",
                 value.name = "F1")

means = apply(df, 2, mean)
means["fold"] = NA
df_means = data.frame("Mean" = means, "Encoding" = names(means))

#svg("overview.svg")
ggplot(df_means, aes(Encoding, Mean)) +
  geom_point(color = "#4C78A8") +
  geom_hline(yintercept=0.8, linetype="dashed", color = "grey") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
#dev.off()

df_filtered = df_melted[df_melted[, "Encoding"] %in% names(means[means>=0.8]), ]
df_filtered$Encoding = as.character(df_filtered[ ,"Encoding"])
df_filtered$Encoding = factor(df_filtered$Encoding,
                              levels = sort(unique(df_filtered$Encoding),
                                            decreasing = TRUE))

#svg("top_encodings.svg")
ggplot(df_filtered, aes(F1, Encoding)) + geom_boxplot(fill = "#4C78A8", color = "black") + theme_bw()
#dev.off()

df_filtered_1 = df_means[df_means[, "Encoding"] %in% names(means[means>=0.5 & means<0.6]), ]
df_filtered_1["group"] = "0.5<=means<0.6"
df_filtered_2 = df_means[df_means[, "Encoding"] %in% names(means[means>=0.6 & means<0.7]), ]
df_filtered_2["group"] = "0.6<=means<0.7"
df_filtered_3 = df_means[df_means[, "Encoding"] %in% names(means[means>=0.7 & means<0.8]), ]
df_filtered_3["group"] = "0.7<=means<0.8"
df_filtered_4 = df_means[df_means[, "Encoding"] %in% names(means[means>=0.8 & means<0.9]), ]
df_filtered_4["group"] = "0.8<=means<0.9"
#df_filtered_5 = df_means[df_means[, "Encoding"] %in% names(means[means>=0.9]), ]
#df_filtered_5["group"] = "means>=0.9"

df_res = rbind(df_filtered_1, df_filtered_2, df_filtered_3, df_filtered_4)#, df_filtered_5)

#svg("violin.svg")
ggplot(df_res, aes(group, Mean)) + geom_violin(trim = FALSE, fill = "#4C78A8") + geom_jitter(height = 0, width = 0.05, size=0.5) + theme_bw()
#dev.off()
