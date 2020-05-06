library(pheatmap)
library(stringr)

setwd("D:/These/NetworkInference/")
load("Data/normalized.count_At.RData")



draw_heatmap <- function(normalized.count, subset = "random", show_rownames = F, title = "Random preview heatmap of log(expression)"){
  
  if (subset == "random") sample_subset <- sample(rownames(normalized.count), size = 100)
  else sample_subset <- subset
  
  mat <- log(normalized.count[sample_subset,] + 1)
  
  sample <- str_split_fixed(colnames(mat), '_',2) [,1]
  samples <- data.frame(sample, row.names = colnames(mat))
  pheatmap(mat, annotation_col = samples, show_rownames = show_rownames, 
           main = title, fontsize = 17)
  
}






draw_distributions <- function(normalized.count, boxplot = T){
  d <- melt(log(normalized.count[sample(rownames(normalized.count), replace = F, size = round(dim(normalized.count)[1]/4, 0)),]+1))
  colnames(d) <- c("gene", "sample", "logCount")
  d$condition <- str_split_fixed(d$sample, "_", 2)[,1]
  
  g <- ggplot(data = d, aes(x=sample, y=log_count))
  
  if (boxplot) {
    g <- g + geom_boxplot(alpha=0.5, lwd=1, aes( fill = condition), outlier.color = "black",outlier.alpha =0.1)
  }else{
    g <- g + geom_violin(alpha=0.5, lwd=1, aes( fill = condition))
  }
  
  g <- g +theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 20),legend.position="bottom",
                legend.title = element_text(size = 20, face="bold"), legend.text = element_text(size=22, angle=0),
                axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 15, angle = -50, hjust = 0, colour = "grey50"),legend.text.align=1,
                axis.title=element_text(size=24)) 

  g
}



normalized.count
