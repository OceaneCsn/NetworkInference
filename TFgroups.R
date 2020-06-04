
load('D:/These/NetworkInference/Data/normalized.count_At.RData')
normalized.count <- normalized.count[,grepl('F', colnames(normalized.count))]


library(DIANE)
library(stringr)

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
top_pos <- pairs[pairs$cor > 0.9,]
top_neg <- pairs[pairs$cor < -0.9,]


net_un <- graph_from_data_frame(top_pos, directed = FALSE)
louvain <- cluster_louvain(net_un)
groups <- membership(louvain)

groups

group_correlated_TFs <- function(normalized.count, regressors, corr_thr = 0.9){
  pairs <- data.frame(t(combn(regressors, 2)))
  pairs$cor <- sapply(paste(pairs[,1], pairs[,2]), get_cor)
  top <- pairs[pairs$cor > corr_thr,]
  
  net_un <- graph_from_data_frame(top, directed = FALSE)
  louvain <- cluster_louvain(net_un)
  groups <- membership(louvain)
  
  new_reg <- c()
  grouped_regs <- data.frame(matrix(nrow = length(unique(groups)), ncol = length(colnames(normalized.count))))
  colnames(grouped_regs) <- colnames(normalized.count)
  rownames(grouped_regs) <- unique(groups)
  
  for(group in unique(groups)){
    tfs <- names(groups)[groups==group]
    mean_tf <- colMeans(normalized.count[tfs,])
    grouped_regs[group,] <- mean_tf
    new_reg <- c(new_reg, paste0("mean_",paste(tfs, collapse = "-")))
  }
  rownames(grouped_regs) <- new_reg
  normalized.count <- rbind.data.frame(normalized.count, grouped_regs)
  return(normalized.count)
}

d <- group_correlated_TFs(normalized.count, regressors)
dim(d$counts)
dim(normalized.count)
co <- d$counts
rownames(co)[rownames(co) %in% regressors]
d$new_regressors

dim(normalized.count)

demo_data_At$gene_info$label[match(tfs, rownames(demo_data_At$gene_info))]




tail(normalized.count)


net <- graph_from_data_frame(top, directed = TRUE)

library(visNetwork)

data <- toVisNetworkData(net_un)

data$edges$value <- data$edges$cor
data$edges$label <- round(data$edges$cor,3)
data$nodes$group <- groups[match(data$nodes$id, names(groups))]
data$nodes$label <- demo_data_At$gene_info$label[match(data$nodes$id, rownames(demo_data_At$gene_info))]

visNetwork(nodes = data$nodes, edges = data$edges) %>% visNodes(font = list("size" = 35))



