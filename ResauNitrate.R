setwd("D:/These/NetworkInference")

load("./Data/DEGsListsFiltered.RData")
load("Data/PlnTFBDRegulatorsList.RData")
load("./Data/normalized.count_At.RData")
load("Data/OntologyAllGenes.RData")

library(GENIE3)
library(bc3net)
library(stringr)
source("Funtions/Network_functions.R")

genes <- DEGs[["cnf Cnf"]]
rownames(ontologies) <- ontologies$ensembl_gene_id

removeNitrateStarv = F; removeIronStarv = F; removeCO2Stress=F;


if (removeIronStarv) normalized.count <- normalized.count[,grepl("F", colnames(normalized.count))]
if (removeNitrateStarv) normalized.count <- normalized.count[,grepl("N", colnames(normalized.count))]
if (removeCO2Stress) normalized.count <- normalized.count[,grepl("c", colnames(normalized.count))]

head(normalized.count)

mat <- genie(normalized.count, regressors = intersect(genes, TF$AGI), targets = genes)
  
  

net = getNet(mat, fixedLinkNumber = 3000)
data <- networkData(net, ontologies, TF)

plotNetwork(data)

save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_LowNitrate_CO2-N-Fe.RData")


plotStats(net)
netStats(net)

data <- dataFer
BC3net <- bc3net(dataset = normalized.count[genes,], alpha1=0.01, alpha2=0.01 )
gcc <- getgcc(BC3net)
data <- networkData(gcc, ontologies, TF)
plotNetwork(data)

save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_FaibleNitrate_CO2-N_BC3NET.RData")



plotStats(net)
data$nodes


############### Reseau nitrate par gaudinier 2018

gaudinier <- read.csv("./Data/GaudinierTF-TargetLinksNitrate.csv", h=T)
colnames(gaudinier)[4] <- "to"
colnames(gaudinier)[1] <- "from"

net <- igraph::graph_from_data_frame(d = gaudinier[,c(1,4)])

E(net)

data <- networkData(net, ontologies)


data$nodes
data$edges
plot.igraph(net,vertex.size=4, vertex.label.cex=0.5)

gaudinier <- data
save(gaudinier, file = "D:/These/NetworkShiny/NetworkData/Gaudinier_Nature_2018_TFs_Nitrates_Ref.RData")


load("D:/These/NetworkShiny/NetworkData/NitrateDEGenes_FortCO2_CO2-N.RData")


data$edges
gaudinier$edges

data$edges$pairs <- paste(data$edges$from, data$edges$to)
gaudinier$edges$pairs <- paste(gaudinier$edges$from, gaudinier$edges$to)
data$edges[intersect(c(data$edges$pairs, "AT4G29080 AT2G14210"), gaudinier$edges$pairs),]

plotNetwork(data)


################### essai avec PLNNetwork

net <- PLN_network(normalized.count, genes)
netStats(net)
save(net, file = "networkPLN_NitrateGenes_NFeCO2.RData")


######################  reseau litterature correlations de gaudinier et Al :


gaudinier_corr <- read.csv("Data/Gaudinier_2018_correlations.csv", h=T)[,c("TF_AGI", "Promoter_AGI")]
net <- igraph::graph_from_data_frame(d = gaudinier_corr)
data <- networkData(net, ontologies)
plotNetwork(data)
save(gaudinier_corr, file = "D:/These/NetworkShiny/NetworkReference/Gaudinier_Correlations_Litterature.RData")
gaudinier_corr <- data
