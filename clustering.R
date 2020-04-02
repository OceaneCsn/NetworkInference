library(igraph)
library("visNetwork")
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
library(plotly)

load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")
load("./Data/OntologyAllGenes.RData")
setwd("D:/These/NetworkInference")
plotNetwork <- function(data){
  visNetwork(nodes = data$nodes, edges = data$edges)%>% 
    #visEdges(smooth = FALSE, color = '#333366') %>% 
    visEdges(smooth = FALSE, arrows = 'to', color = '#333366') %>% 
    visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
               maxVelocity = 10, stabilization = F)%>%
    visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = TRUE)%>% 
    visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
    visGroups(groupname = "Regulator", size = 28,
              color = list("background" = "#003399", "border"="#FFFFCC"), shape = "square") %>% 
    visGroups(groupname = "Target Gene", color = "#77EEAA") %>% 
    visNodes(borderWidth=0.5, font=list("size"=36)) 
}


net <- igraph::graph_from_data_frame(d = data$edges, directed = F)
communities <- cluster_louvain(net)
communities$membership
membership <- membership(communities)
data$nodes$group <- membership[match(data$nodes$id, names(membership))]
plotNetwork(data)

data$edges


universe <- ontologies$entrezgene_id


OntologyEnrich <- function(ids, universe, plot = T, simCutoff = 0.8){
  # ids and universe must be entrez gene ids
  ego <- enrichGO(gene = ids,
                  OrgDb = org.At.tair.db,
                  ont = "BP",
                  universe = universe,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)
  
  # Elimine les redondances, en fusionnant les GO terms dont la similarite excede le cutoff
  simpOnt <- clusterProfiler::simplify(ego, cutoff=simCutoff, by="p.adjust", select_fun=min)
  result <- simpOnt@result
  if(plot){
    print(barplot(simpOnt, showCategory = 40, font.size = 10))
    print(emapplot(simpOnt, layout = "kk"))
  }
  return(simpOnt)
}

cluster <- membership[1]
ids <- as.character(ontologies[match(data$nodes[data$nodes$group==cluster,]$id, ontologies$ensembl_gene_id),]$entrezgene_id)
simOnt <- OntologyEnrich(ids, as.character(universe))
simOnt@result <- simOnt@result[order(-simOnt@result$p.adjust),]
simOnt@result$GeneRatio
ggplotly(ggplot(data= simOnt@result, aes(x=Description, y=GeneRatio, fill=p.adjust)) + geom_bar(stat = "identity")+ coord_flip() +
           ggtitle(paste("Enriched Ontologies for module", cluster)) + ylab("") + xlab("Gene Ratio")+
           theme(plot.title = element_text(size=15, face="bold"),
                 legend.title = element_text(size = 20, face="bold"), legend.text = element_text(size=15),
                 axis.text.y = element_text(size = 15, angle = 20), axis.text.x = element_text(size =8, angle = 0, hjust = 0, colour = "grey50"),
                 axis.title.y=element_blank()) )


idsList <- list()
for(k in unique(membership)){
  idsList[[as.character(k)]] <- na.omit(as.character(ontologies[match(data$nodes[data$nodes$group==k,]$id, ontologies$ensembl_gene_id),]$entrezgene_id))
}

compareOnt <- function(idsList, universe, simCutoff = 0.8){
  ckreg <- compareCluster(geneCluster = idsList, fun = "enrichGO", OrgDb = org.At.tair.db, ont = "BP", pAdjustMethod = "BH", 
                          pvalueCutoff = 0.01, qvalueCutoff = 0.05, universe = as.character(universe))
  ckreg@compareClusterResult
  simCk <- clusterProfiler::simplify(ckreg, cutoff=simCutoff, by="p.adjust", select_fun=min)
  resCk <- simCk@compareClusterResult
  print(clusterProfiler::dotplot(simCk, x = ~Cluster, showCategory = 25, font.size = 7))
  return(resCk)
}
compareOnt(idsList=idsList, universe)



