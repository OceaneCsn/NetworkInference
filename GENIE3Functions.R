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
removeNitrateStarv = F; removeIronStarv = T

if (removeIronStarv) normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
if (removeNitrateStarv) normalized.count <- normalized.count[,grepl("N", colnames(normalized.count))]



#' Inserts the random variables
#'
#' @param normalized.count expression data to which we wan't to add the spy variables
#' @param regressors regulators to use to generate the spies
#' @param n_spies number of random varaibles
#' @param keepColMeans weather or not to use means specific to the conditions or a global mean
#'
#' @return the expression dataframe with spy variables inside
#' @export
insertSpyVariable <- function(normalized.count, regressors = row.names(normalized.count), 
                              n_spies = 3, keepColMeans = T){
  means <- colMeans(normalized.count[regressors,])
  
  spies <- data.frame(normalized.count[sample(rownames(normalized.count), size=n_spies),], 
                      row.names = paste0("spy_",seq(1:n_spies)))
  for(cond in colnames(normalized.count)){
    if(keepColMeans) spies[,cond] <- rpois(lambda = means[cond], n = n_spies)
    else spies[,cond] <- rpois(lambda = mean(means), n = n_spies)
  }
  normalized.count <- rbind.data.frame(normalized.count, spies)
  print(tail(normalized.count))
  return(normalized.count)
}

insertSpyVariable(normalized.count)
#s2 <- insertSpyVariable(normalized.count, keepColMeans = F)

#' Function to combine genie3 and our thresholding method
#'
#' @param normalized.count expression data
#' @param regressors regulators (TFs) names in normalized.count
#' @param targets target genes
#' @param nTrees number of trees in GENIE3
#' @param nCores number of cores in GENIE3
#' @param pval quantile to use on the null distribution
#' @param varNorm weather to normalize the expression data to have 
#' a unit variance
#' @param n_spies number of random variables
#' @param keepColMeans weather or not to use means specific to the conditions or a global mean
#' @param returnNlinks weather or not to return only the number of links, not the network
#'
#' @return the infered graph, or the number of links of the inferred graph
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


genie3_gene_specific <- function(normalized.count, regressors=NA, targets=NA, nTrees=5000, nCores=5, 
                   pval = 0.05, varNorm = F, n_spies = 5, keepColMeans = T, returnNlinks = F){
  
  normalized.count <- insertSpyVariable(normalized.count, regressors,
                                        n_spies = n_spies, keepColMeans = keepColMeans)
  
  regressors <- intersect(rownames(normalized.count),regressors)
  regressors = c(regressors, rownames(normalized.count)[grepl('spy_', rownames(normalized.count))])
  
  
  if (varNorm) normalized.count <- normalized.count/rowSds(normalized.count)
  
  
  
  mat <- GENIE3(as.matrix(normalized.count), regulators = regressors, targets = targets, 
                treeMethod = "RF", K = "sqrt", nTrees = nTrees, nCores = nCores, verbose = T)
  
  thrs <- mat[grepl('spy_', rownames(mat)),]
  thrs <- t(thrs)
  maxs <- apply(thrs,1,max)

  for(target in rownames(mat)){
    mat[target,] <- ifelse(mat[target,] > maxs[target], mat[target,],-1)
  }
  
  links <- getLinkList(mat, threshold = 0)
  print(dim(links)[1])
  if (returnNlinks) return(dim(links)[1]) 
  print(paste0("Number of links : ", dim(links)[1]))
  g <- graph_from_data_frame(links, directed = T)
  return(g)
  
}


net <- genie3_gene_specific(normalized.count,regressors = intersect(TF$AGI, genes), targets = genes, varNorm = F, pval = 0.001,
              n_spies = 30, keepColMeans = F)


maxs <- net

head(maxs)
df <- data.frame(genes = names(maxs), maxImp = maxs)
df$Sd <- rowSds(normalized.count[match(df$genes, rownames(normalized.count)),])
df$mean <- rowMeans(normalized.count[match(df$genes, rownames(normalized.count)),])

ggplot(data = df, aes(y = maxImp, x = Sd/mean)) + geom_hex() + ggtitle("importance max de la variable aléatoire en fonction du ratio sd/mean")

cor(df$maxImp, df$Sd/df$mean, method = "spearman")
#10 % de correlation entre l'importance et le ratio sd mean

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
          #print(meanSpecific)
          nlinks <- genie3(normalized.count,regressors = intersect(TF$AGI, genes), targets = genes, varNorm = F, pval = 0.001,
                                 n_spies = n_spies, keepColMeans =meanSpecific, returnNlinks = T)
          res[cpt, ] = c(meanSpecific, norm_var, pvalue, n_spies, rep, nlinks)
          print("res here :")
          print(c(meanSpecific, norm_var, pvalue, n_spies, rep, nlinks))
          cpt = cpt + 1
        }
      }
    }
  }
}

library(ggplot2)


res$meanSpecific <- ifelse(res$meanSpecific==1,"condition specific mean", "global mean")
res$norm_var <- ifelse(res$norm_var==1, "normalized variance", "not normalized variance")

ggplot(data = res, aes(x = n_spies , y = n_links, color = as.factor(pvalue)))+facet_wrap(c("norm_var", "meanSpecific")) + geom_point(size = 4, alpha =0.6) + ylim(c(0,5000))

save(res, file = "D:/These/NetworkInference/res_exploparam_spyVariables.RData")

############# exploration paramétrique gene specific

n <- 6*2*4
res <- data.frame(meanSpecific = seq(1:n), 
                  n_spies = seq(1:n),rep = seq(1:n), n_links = seq(1:n))
cpt = 1
for (meanSpecific in c(TRUE, FALSE)){
    for(n_spies in c(3,10,20,40)){
      for(rep in seq(1:6)){
        print(meanSpecific)
        nlinks <- genie3_gene_specific(normalized.count,regressors = intersect(TF$AGI, genes), targets = genes, varNorm = F, pval = 0.001,
                         n_spies = n_spies, keepColMeans =meanSpecific, returnNlinks = T)
        res[cpt, ] = c(meanSpecific, n_spies, rep, nlinks)
        print(c(meanSpecific, n_spies, rep, nlinks))
        cpt = cpt + 1
      }
    }
    
}


library(ggplot2)


res$meanSpecific <- ifelse(res$meanSpecific==1,"condition specific mean", "global mean")
res$norm_var <- ifelse(res$norm_var==1, "normalized variance", "not normalized variance")
min(res$n_links)
# au minimum on a 21200 noeuds 

ggplot(data = res, aes(x = n_spies , y = n_links, color = n_spies))+facet_wrap("meanSpecific") + geom_point(size = 4, alpha =0.6) 

save(res, file = "D:/These/NetworkInference/res_exploparam_spyVariables.RData")





# etude des correlations nb cibles et correlation de paires de TFs
load("D:/These/NetworkInference/corrTFsCommonTargets.RData")

top$commonRatio <- top$CommonTargetsNumber/apply(top[,c("tf1_target_Number", "tf2_target_Number")], 1, min)
library(ggplot2)
ggplot(data = top, aes(y = commonRatio, x = abs.cors.)) + geom_hex()


###################################### genie3 sur les cibles des deux TFs les plus corrélés WRKY75 et NAC055

getTargetsAgi <- function(tf, data){
  targets <- data$edges[data$edges$from==tf,"to"]
  names <- data$nodes[data$nodes$id %in% targets, "id"]
  return(names)
}


load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")
targets <- union(getTargetsAgi(tf = "AT5G13080", data), getTargetsAgi(tf = "AT3G15500", data))
# on a 62 cibles au total, reparties entre les deux TFs

########## en mettant les deux TFs

targets_WRKY <- getTargetsAgi(tf = "AT5G13080", data)
targets_NAC <- getTargetsAgi(tf = "AT3G15500", data)
common <- intersect(targets_WRKY, targets_NAC)


######### en enlevant WRKY75, on refait un reseau avec 3000
library(visNetwork)
library(stringr)
# net was generated in genie3.rmd sale je sais
data_sans_NAC55 <- networkToData(net, ontologies, TF)
targets_NAC_alone <- getTargetsAgi(tf = "AT3G15500", data_sans_WRKY)

# NAC passe de 36 à 43 cibles. Voyons le recoupement avec celles de WRK dans le reseau classique : 
length(intersect(targets_NAC_alone, targets_WRKY))
# parmis les 43 cibles de NAC, 13 etaient aussi des cibles de WRKY dans le reseau classique (contre 9 de partagées dans le reseau classique)
# on peut donc supposer que NAC a volé 4 cibles a WRKY en son abscence (mais largement pas toutes, il en reste 22 autres)
# si on regarde ces autres cibles : 

common_sans_WRKY <- intersect(targets_NAC_alone, targets_WRKY)
not_stolen_by_NAC <- targets_WRKY[!targets_WRKY %in% common_sans_WRKY]


getRegulators_agi <- function(genes, data){
  regs <- data$edges[data$edges$to %in% genes,"from"]
  names <- data$nodes[data$nodes$id %in% regs, "id"]
  return(unique(names))
}

# on va regarder la correlation des TFs qui ont volé des cibles a WRK avec WRKY
getRegulators_agi(not_stolen_by_NAC, data_sans_WRKY)


getCor <- function(tf){
  return(abs(cor(normalized.count[tf,], normalized.count["AT5G13080",], method = "spearman")))
}
cors <- sapply(getRegulators_agi(not_stolen_by_NAC, data_sans_WRKY), getCor)
hist(cors, breaks = 20)
summary(cors)



######### en mettant NAC055

targets_WRKY_alone <- getTargetsAgi(tf = "AT5G13080", data_sans_NAC55)
# c'est fort etonnant!! quand Nac n'est pas la, WRKY fait moins de connexions, on passe de 35 à 30...

length(intersect(targets_WRKY_alone, targets_WRKY))
# on a perdu 3 connexions que faisait WRKY


length(intersect(targets_WRKY_alone, targets_NAC))
# 8 de ces cibles font partie de celles de NAC avait dans le reseau classique

common_sans_NAC <- intersect(targets_WRKY_alone, targets_NAC)
not_stolen_by_WRKY <- targets_NAC[!targets_NAC %in% common_sans_NAC]
#les cibles de NAC que WRKY n'a pas recuperees

cors <- sapply(getRegulators_agi(not_stolen_by_WRKY, data_sans_NAC55), getCor)
hist(cors, breaks = 20)
summary(cors)

# pour ceux qui récupèrent des cibles de NAC, on n'a meme pas de forte correlation avec NAC, ça pourrait etre que les combinatoires impliquent plusieurs facteurs