
load('D:/These/NetworkInference/Data/normalized.count_At.RData')
normalized.count <- normalized.count[,grepl('F', colnames(normalized.count))]


library(DIANE)
data("demo_data_At")


pairs <- data.frame(t(combn(regressors, 2)))

get_cor <- function(pair){
  tf1 <- str_split_fixed(pair, ' ',2)[1]
  tf2 <- str_split_fixed(pair, ' ',2)[2]
  return(cor(normalized.count[tf1,], normalized.count[tf2,], method = "spearman"))
}

pairs$cor <- sapply(paste(pairs[,1], pairs[,2]), get_cor)
head(pairs)


top <- pairs[abs(pairs$cor) > 0.9,]


net <- graph_from_data_frame(top, directed = TRUE)

library(visNetwork)
library(stringr)

data <- toVisNetworkData(net)

data$edges$value <- data$edges$cor
data$edges$label <- round(data$edges$cor,3)

visNetwork(nodes = data$nodes, edges = data$edges)

net_un <- graph_from_data_frame(top, directed = FALSE)
louvain <- cluster_louvain(net_un)
groups <- membership(louvain)


data$nodes$group <- groups[match(data$nodes$id, names(groups))]
visNetwork(nodes = data$nodes, edges = data$edges)




