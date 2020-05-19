

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
data$nodes


#ajout du cluster louvain
net <- igraph::graph_from_data_frame(d = data$edges, directed = F)
communities <- cluster_louvain(net)
communities$membership
membership <- membership(communities)
data$nodes$louvainCluster <- membership[match(data$nodes$id, names(membership))]

newData <- data
newData$nodes$group <- newData$nodes$coseqCluster
plotNetwork(newData)


# ajout du cluster coseq 
load("D:/These/ClusteringAnalysis/Clusterings/AmbientCO2_LowNitrateFe-ElevatedCO2_LowNitrateFeNoIronStarv.RData")

#load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_LowNitrate_CO2-N .RData")
data$nodes$coseqCluster <- cluster[[1]][match(data$nodes$id, names(cluster[[1]]))]



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
save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N_PvalueRandomVariable0.1.RData" )

colnames(data$nodes)
data$nodes

data$nodes$group <- data$nodes$consensusCluster
data$edges$color <- "grey"
plotNetwork(data)%>%  visGroups(groupname = "-1", color = "lightgrey")
save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData" )



############################ corrélation des TFs entre eux
load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")

load("D:/These/NetworkShiny/Data/normalized.count_At.RData")

edges <- paste(data$edges$from, data$edges$to)
reversedEdges <- paste(data$edges$to, data$edges$from)

intersect(edges, reversedEdges)

cors <- c()
names <- c()
#corEdgesDouble <- c()
tfs <- data$nodes[data$nodes$group == "Regulator","id"]
for (tf1 in tfs) {
  for (tf2 in tfs) {
    if (tf1 != tf2) {
      cors <- c(cors, cor(normalized.count[tf1, ], normalized.count[tf2, ], method = "spearman"))
      names <- c(names, paste(tf1, tf2))
      # if(paste(tf1, tf2) %in% recip){
      #   corEdgesDouble <- c(corEdgesDouble, cor(normalized.count[tf1, ], normalized.count[tf2, ], method = "spearman"))
      # }
    }
  }
}

d <- data.frame(names, abs(cors))
d <- d[seq(1,dim(d)[1],2),]

#d <- d[!duplicated(d$abs.cors.),]


head(d[order(-d$abs.cors.),])

# si on prend les top 0.5% des correlations pairwises entre tfs : 
q <- quantile(d$abs.cors., probs = 0)
top <- d[d$abs.cors.> q,]
top$tf1 <- str_split_fixed(top$names, " ", 2)[,1]
top$tf2 <- str_split_fixed(top$names, " ", 2)[,2]


load("D:/These/NetworkShiny/Data/OntologyAllGenes.RData")

top$tf1_name <- ontologies[match(top$tf1, ontologies$ensembl_gene_id),"external_gene_name"]
top$t21_name <- ontologies[match(top$tf2, ontologies$ensembl_gene_id),"external_gene_name"]


getTargets <- function(tf){
  targets <- data$edges[data$edges$from==tf,"to"]
  names <- data$nodes[data$nodes$id %in% targets, "label"]
  #return(paste(names, collapse = ', '))
  return(names)
}

getCommonTargets <- function(pair){
  tf1 <- str_split_fixed(pair, ' ', 2)[,1]
  tf2 <- str_split_fixed(pair, ' ', 2)[,2]
  
  # print(getTargets(tf1))
  # print(getTargets(tf2))
  # 
  print(intersect(getTargets(tf1), getTargets(tf2)))
  
  return(length(intersect(getTargets(tf1), getTargets(tf2))))
}

getTargetsNumber <- function(tf){
  return(length(data$edges[data$edges$from==tf,"to"]))
}

top$tf1_target_Number <- sapply(top$tf1, getTargetsNumber)
top$tf2_target_Number <- sapply(top$tf2, getTargetsNumber)

#top$tf1_targets <- sapply(top$tf1, getTargets)
#top$tf2_targets <- sapply(top$tf2, getTargets)


tf1 = "AT3G25730"
tf2 = "AT1G77450"


top$CommonTargetsNumber <- sapply(top$names, getCommonTargets)

#hist(cors, breaks = 30)
#hist(corEdgesDouble, breaks = 30)

plot(x = top$abs.cors., y = top$CommonTargetsNumber, main = "Nombre de cibles communes en fonction de la correlation pour toutes les paires de TFs")
cor(top$abs.cors., top$CommonTargetsNumber, method = "spearman")

save(top, file = "D:/These/NetworkInference/corrTFsCommonTargets.RData")
