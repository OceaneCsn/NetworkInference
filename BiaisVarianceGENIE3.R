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





load('D:/These/NetworkInference/corrTFsCommonTargets.RData')
pair <- c("AT5G13080", "AT3G15500")
tf1 <- pair[1]
tf2 <- pair[2]

source("Funtions/Network_functions.R")

load("./Data/DEGsListsFiltered.RData")
load("./Data/PlnTFBDRegulatorsList.RData")
load("./Data/normalized.count_At.RData")
load("./Data/OntologyAllGenes.RData")


pval = 0.001
# DE genes
genes <- DEGs[["cnF CnF"]]
print(length(genes))

# expression data
normalized.count <- normalized.count[,grepl('F', colnames(normalized.count))]
head(normalized.count)

# TFs
regressors = intersect(TF$AGI, genes)



mat <- GENIE3(normalized.count, targets = genes, regulators = regressors)
links <- getLinkList(mat, reportMax = 3000)



load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")

net <-  graph_from_data_frame(data$edges)

tf1_targets <- neighbors(net, tf1, mode = "out")
tf2_targets <- neighbors(net, tf2, mode = "out")

common <- intersect(tf1_targets, tf2_targets)

c1 <- common[4]

tf1_only <- tf1_targets[!tf1_targets %in% common]
tf2_only <- tf2_targets[!tf2_targets %in% common]


plot(d[,1], d[,2], main = "Importance des genes cibles communs pour tf1 et tf2")

# inferences sans l'un des deux

e <- data$edges

mat_sans_1 <- GENIE3(normalized.count, targets = genes, regulators = regressors[regressors!=tf1])
mat_sans_2 <- GENIE3(normalized.count, targets = genes, regulators = regressors[regressors!=tf2])


# resultats

d <- data.frame("importance TF1" = mat[tf1, common], "importance TF2" = mat[tf2, common], 
                "importance TF1 sans TF2" = mat_sans_2[tf1, common],
                "importance TF2 sans TF1" = mat_sans_1[tf2, common])
d$ratioTF1 <- d[,3]/d[,1]
d$ratioTF2 <- d[,4]/d[,2]
d$target_type <- "common target" 


d_1 <- data.frame("importance TF1" = mat[tf1, tf1_only], "importance TF2" = mat[tf2, tf1_only], 
                  "importance TF1 sans TF2" = mat_sans_2[tf1, tf1_only],
                  "importance TF2 sans TF1" = mat_sans_1[tf2, tf1_only])
d_1$ratioTF1 <- d_1[,3]/d_1[,1]
d_1$ratioTF2 <- d_1[,4]/d_1[,2]
d_1$target_type <- "TF1 target only" 


d_2 <- data.frame("importance TF1" = mat[tf1, tf2_only], "importance TF2" = mat[tf2, tf2_only], 
                  "importance TF1 sans TF2" = mat_sans_2[tf1, tf2_only],
                  "importance TF2 sans TF1" = mat_sans_1[tf2, tf2_only])
d_2$ratioTF1 <- d_2[,3]/d_2[,1]
d_2$ratioTF2 <- d_2[,4]/d_2[,2]
d_2$target_type <- "TF2 target only" 


res <- rbind.data.frame(d, d_1, d_2)


library(ggplot2)


ggplot(data = res, aes(x = target_type, y = ratioTF1, fill = target_type)) + 
  geom_boxplot(alpha = 0.3) + ggtitle("Ratio entre l'importance du TF1 sans TF2, et avec TF2")

ggplot(data = res, aes(x = target_type, y = ratioTF2, fill = target_type)) + 
  geom_boxplot(alpha = 0.3) + ggtitle("Ratio entre l'importance du TF2 sans TF1, et avec TF1")

ggplot(data = res, aes(x = target_type, y = importance.TF1/importance.TF2, fill = target_type)) + 
  geom_boxplot(alpha = 0.3) + ggtitle("Ratio entre l'importance du TF1 et TF2 dans les données complètes")



