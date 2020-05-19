########## test biais de variance pour inference GENIE3 

load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")

load("D:/These/NetworkShiny/Data/normalized.count_At.RData")

net <- igraph::graph_from_data_frame(d = data$edges)
library(igraph)

in_degrees <- igraph::degree(net, mode = "in")
d <- data.frame(degree = as.numeric(in_degrees), genes = names(in_degrees))

normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
head(normalized.count)

getVar <- function(gene){
  return(var(normalized.count[gene,]))
}
d$var = as.numeric(sapply(d$genes, getVar))
d <- d[order(-d$var),]

library("ggplot2")
#ggplot(data = d, aes(x= var, y=degree)) + geom_hex(bins = 20) +xlim(c(0,52029485.41))
ggplot(data = d, aes(x= sqrt(var), y=degree)) + geom_hex() + xlim(c(0,7500)) +ggtitle("Gene in-degree depending on gene standard deviation")

# pearson
cor(d$degree, d$var)
#spearman
cor(d$degree, d$var, method = "spearman")


out_degrees <- igraph::degree(net, mode = "out", v = d$gene)
d$out_degree <- out_degrees[match(d$genes, names(out_degrees))]
ggplot(data = d, aes(x= sqrt(var), y=out_degree)) + geom_point() + xlim(c(0,7500)) +ggtitle("Gene out-degree depending on gene standard deviation")
summary(lm(d$degree~d$var))


## variance des genes qui n'ont pas ete pris pour l'inference

load("Data/DEGsListsFiltered.RData")
length(DEGs$`cnF CnF`)
length(data$nodes$id)

df <- data.frame(genes = DEGs$`cnF CnF`)
df$in_network <- ifelse(df$genes %in% data$nodes$id, 1, 0)
df$var <- as.numeric(sapply(df$genes, getVar))


ggplot(data = df, aes(x = as.factor(in_network), y = var, fill = as.factor(in_network))) + geom_boxplot(alpha = 0.4) + ylim(0,1000) +
  stat_compare_means(method = "t.test")






