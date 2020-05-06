library(stringr)
library(GENIE3)
library(igraph)
library(visNetwork)
library(genefilter)

setwd("D:/These/NetworkInference")

load("./Data/DEGsListsFiltered.RData")
load("Data/PlnTFBDRegulatorsList.RData")
load("./Data/normalized.count_At.RData")
load("Data/OntologyAllGenes.RData")


source("Funtions/Network_functions.R")

genes <- DEGs[["cnF CnF"]]
#rownames(ontologies) <- ontologies$ensembl_gene_id
removeNitrateStarv = F; removeIronStarv = T

if (removeIronStarv) normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
if (removeNitrateStarv) normalized.count <- normalized.count[,grepl("N", colnames(normalized.count))]




insertSpyVariable <- function(normalized.count, regressors = row.names(normalized.count), 
                              n_spies = 5, keepColMeans = T){
  means <- colMeans(normalized.count[regressors,])
  
  spies <- data.frame(normalized.count[sample(rownames(normalized.count), size=n_spies),], row.names = paste0("spy_",seq(1:n_spies)))
  for(cond in colnames(normalized.count)){
    if(keepColMeans) spies[,cond] <- rpois(lambda = means[cond], n = n_spies)
    else spies[,cond] <- rpois(lambda = mean(means), n = n_spies)
  }
  
  normalized.count <- rbind.data.frame(normalized.count, spies)
  print(tail(normalized.count))
  return(normalized.count)
}

#s <- insertSpyVariable(normalized.count)
#s2 <- insertSpyVariable(normalized.count, keepColMeans = F)

genie3 <- function(normalized.count, regressors=NA, targets=NA, nTrees=5000, nCores=5, 
                  pval = 0.05, varNorm = F, n_spies = 5, keepColMeans = T, returnNlinks = F){
  
  normalized.count <- insertSpyVariable(normalized.count, regressors,
                                        n_spies = n_spies, keepColMeans = keepColMeans)
  
  regressors <- intersect(rownames(normalized.count),regressors)
  regressors = c(regressors, rownames(normalized.count)[grepl('spy_', rownames(normalized.count))])

  
  if (varNorm) normalized.count <- normalized.count/rowSds(normalized.count)
  

  
  mat <- GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores, verbose = T)
  
  thrs <- c(mat[grepl('spy_', rownames(mat)),])
  threshold = quantile(thrs, probs = 1-pval)
  
  
  print(threshold)

  
  links <- getLinkList(mat, threshold = threshold)
  links <- links[!grepl('spy_', links[,1]),]
  if (returnNlinks) return(dim(links)[1]) 
  print(paste0("Number of links : ", dim(links)[1]))
  g <- graph_from_data_frame(links, directed = T)
  return(g)
}


net <- genie3(normalized.count,regressors = intersect(TF$AGI, genes), targets = genes, varNorm = F, pval = 0.001,
              n_spies = 3, keepColMeans = F)


networkToData <- function(net, ontologies, TF){
  data <- toVisNetworkData(net)
  data$nodes$Ontology <- ontologies[match(data$nodes$id, ontologies$ensembl_gene_id), "external_gene_name"]
  data$nodes$description <- ontologies[match(data$nodes$id, ontologies$ensembl_gene_id), "description"]
  data$nodes$group <- ifelse(data$nodes$id %in% TF$AGI, "Regulator", "Target Gene")
  
  data$nodes$label <- data$nodes$Ontology
  data$edges$color <- '#333366'
  data$edges$value <- data$edges$weight
  data$edges$Regulator_Name <-
    ontologies[match(data$edges$from, ontologies$ensembl_gene_id), ]$external_gene_name
  data$nodes$description <-
    ifelse(
      grepl("\\[", data$nodes$description),
      str_split_fixed(data$nodes$description, "\\[", 2)[, 1],
      data$nodes$description
    )
  data$edges$Target_Name <-
    ontologies[match(data$edges$to, ontologies$ensembl_gene_id), ]$external_gene_name
  return(data)
}

data <- networkToData(net, ontologies, TF)


save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N_PvalueRandomVariable0.01.RData")

plotNetwork(data)
dim(data$nodes)


dataTF <- function(data){
  newData <- data
  newData$nodes <- newData$nodes[newData$nodes$group == "Regulator", ]
  dim(newData$nodes)
  newData$edges <-
    newData$edges[newData$edges$from %in% newData$nodes$id &
                    newData$edges$to %in% newData$nodes$id , ]
  return(newData)
}
plotNetwork(dataTF(data))



############# erretes doubles

edges <- paste(data$edges$from, data$edges$to)
revEdges <- paste(data$edges$to, data$edges$from)
length(intersect(edges, revEdges))/dim(dataTF(data)$edges)[1]

recip <- intersect(edges, revEdges)


############# exploration paramétrique

n <- 6*2*2*3*3
res <- data.frame(meanSpecific = seq(1:n), norm_var = seq(1:n),pvalue = seq(1:n),
                  n_spies = seq(1:n),rep = seq(1:n), n_links = seq(1:n))
cpt = 1
for (meanSpecific in c(TRUE, FALSE)){
  for (norm_var in c(TRUE, FALSE)){
    for(pvalue in c(0.05,0.01,0.001)){
      for(n_spies in c(3,10,20)){
        for(rep in seq(1:6)){
          print(meanSpecific)
          nlinks <- genie3(normalized.count,regressors = intersect(TF$AGI, genes), targets = genes, varNorm = norm_var, pval = pvalue,
                                 n_spies = n_spies, keepColMeans =meanSpecific, returnNlinks = T)
          res[cpt, ] = c(meanSpecific, norm_var, pvalue, n_spies, rep, nlinks)
          print(c(meanSpecific, norm_var, pvalue, n_spies, rep, nlinks))
          cpt = cpt + 1
        }
      }
    }
  }
}
save(res, file = "D:/These/NetworkInference/res_exploparam_spyVariables.RData")
