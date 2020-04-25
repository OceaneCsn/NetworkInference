

setwd("D:/These/NetworkInference")

library(visNetwork)
library(ggplot2)
source("Funtions/Network_functions.R")


load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_LowNitrate_CO2-N-Fe.RData")
load("./Data/NitrateGenes.RData")

data$nodes

#ajout du score nitrate
res <- data.frame(Gene = data$nodes$id)
for(paper in colnames(nGenes)){
  res[,paper] <- ifelse(res$Gene %in% toupper(nGenes[,paper]), 1, 0)
}
res$NitrateScore <- rowSums(res[,grepl("_", colnames(res))])
res$Gene
data$nodes$NitrateScore <- res$NitrateScore
data$nodes

#ajout du cluster louvain
net <- igraph::graph_from_data_frame(d = data$edges, directed = F)
communities <- cluster_louvain(net)
communities$membership
membership <- membership(communities)
data$nodes$louvainCluster <- membership[match(data$nodes$id, names(membership))]

newData <- data
newData$nodes$group <- newData$nodes$louvainCluster
plotNetwork(newData)

save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_LowNitrate_CO2-N-Fe.RData" )

# Nombre de cibles pour un TF
degree(net)
tf = "AT2G33710"
getTargets <- function(data, tf){
  targets <- data$edges[data$edges$from==tf,"to"]
  names <- data$nodes[data$nodes$id %in% targets, "label"]
  return(paste(names, collapse = ';'))
}

getTargetsNumber <- function(data, tf){
  return(length(data$edges[data$edges$from==tf,"to"]))
}

################ sorting TFs sor cluster k
k = 1
clusterType = "louvainCluster"
topRate = 20
dataK <- data$nodes[data$nodes[,clusterType]==k & data$nodes$group=="Regulator",]

top <- dataK[order(-dataK$ranking),]
top <- top[1:round(dim(dataK)[1]*topRate/100, 0),]
top$targetNumber <- sapply(top$id, getTargetsNumber)
top$targets <- sapply(top$id, getTargets)
top[, c("label", "description", "ranking", "targetNumber", "targets", "NitrateScore")]
#data$nodes$description <- ifelse(grepl("\\[", data$nodes$description), str_split_fixed(data$nodes$description, "\\[", 2)[,1], data$nodes$description)


############ comparer louvain et coseq

library(ggplot2)
library(corrplot)
load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData")


louvain <- data$nodes$louvainCluster
names(louvain) <- data$nodes$id
coseq <- data$nodes$coseqCluster
names(coseq) <- data$nodes$id

mat <- matrix(nrow = max(louvain), ncol=max(coseq))
rownames(mat) <- 1:max(louvain)
colnames(mat) <- 1:max(coseq)
for(louv in 1:max(louvain)){
  louvGenes <- names(louvain[louvain==louv])
  for(cos in 1:max(coseq)){
    cosGenes <- names(coseq[coseq==cos])
    mat[louv,cos] <- length(intersect(louvGenes, cosGenes))
  }
}

heatmap(mat)

library(ramify)
# louvain clusters in names, values is the coseq cluster with max shared genes
maxCoseqForLouvain <- argmax(mat, row=T)

consensus <- function(gene, data){
  coseqClust <- data$nodes[gene,"coseqCluster"]
  louvainClust <- data$nodes[gene,"louvainCluster"]
  if(coseqClust == maxCoseqForLouvain[louvainClust]){
    return(louvainClust)
  }
  else return(-1)
}

data$nodes$consensusCluster <- sapply(data$nodes$id, consensus, data)

data$nodes$group <- data$nodes$consensusCluster
data$edges$color <- "grey"
plotNetwork(data)%>%  visGroups(groupname = "-1", color = "lightgrey")
save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData" )
