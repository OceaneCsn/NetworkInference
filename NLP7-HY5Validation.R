setwd("D:/These/NetworkInference")

library(visNetwork)
library(ggplot2)
source("Funtions/Network_functions.R")


############################################ Cibles NLP7 dans notre reseau


load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData")

query <- dbSendQuery(con, "SELECT * FROM binding WHERE tf_agi = ?")
dbBind(query, list( "AT4G24020"))
res <- data.frame(dbFetch(query))
targets <- res$gene_agi


data$nodes$NLP_Target <- ifelse(data$nodes$id %in% targets, 1, 0)

save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData")

sum(data$nodes$NLP_Target)
sum(data$nodes$NLP_Target)/dim(data$nodes)[1]


data$nodes$group <- data$nodes$NLP_Target

plotNetwork(data)%>%visGroups(groupname = "1", color = "darkred")


load("./Data/OntologyAllGenes.RData")
length(targets)/dim(ontologies)[1]


######################################### Cibles Chip Seq dans notre réseau

load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData")

data$nodes

########## data chip 2007
hy5 <- read.csv("Data/HY5TargetsLee2007.csv", h = T, sep = '\t', stringsAsFactors = F)
hy5Targets <- toupper(hy5$GeneID)
data$nodes$HY5_Target_Chip <- ifelse(data$nodes$id %in% hy5Targets, 1, 0)


hy5edges <- data$edges[data$edges$from == "AT5G11260",]


data$edges$HY5_Chip_Validated <- ifelse(data$edges$to %in% hy5Targets & data$edges$from == "AT5G11260", 1, 0)
sum(data$edges$HY5_Chip_Validated)
18/90

data$nodes$group <- data$nodes$HY5_Target_Chip
sum(data$nodes$HY5_Target_Chip)/dim(data$nodes)[1]

data$edges
plotNetwork(data)


hy5Dap <- hy5edges[hy5edges$DapSeqAproved,"to"]
hy5Chip <- hy5edges[hy5edges$HY5_Chip_Validated==1,"to"]

common <- intersect(hy5Chip, hy5Dap)
commonHy5 <- data$nodes[data$nodes$id %in% common,]

sum(data$edges$to %in% hy5Targets & data$edges$from == "AT5G11260" )

########### data 2018 circadian rythm

hy5BL <- read.csv("Data/Chip-HY5-2018-Circadian-BL.csv", sep = '\t', h=T, stringsAsFactors = F)
hy5BL <- unique(hy5BL$Gene.ID)

hy5RL <- read.csv("Data/Chip-HY5-2018-Circadian-RL.csv", sep = '\t', h=T, stringsAsFactors = F)
hy5RL <- unique(hy5RL$Gene.ID)

#hy5_2018 <- intersect(hy5RL, hy5BL)
hy5_2018 <- union(hy5RL, hy5BL)

############# intersection des deux

#hy5 <- intersect(hy5Targets, hy5_2018)
hy5 <- union(hy5Targets, hy5_2018)


# combien sont confirmées par les papiers
data$edges$HY5_Chip_Validated <- ifelse(data$edges$to %in% hy5 & data$edges$from == "AT5G11260", 1, 0)
sum(data$edges$HY5_Chip_Validated)


# combien de ces interactions sont aussi validées par DapSeq?
hy5edges <- data$edges[data$edges$from == "AT5G11260",]

hy5Dap <- hy5edges[hy5edges$DapSeqAproved,"to"]
hy5Chip <- hy5edges[hy5edges$HY5_Chip_Validated==1,"to"]

common <- intersect(hy5Chip, hy5Dap)
commonHy5 <- data$nodes[data$nodes$id %in% common,]

hist(hy5edges$weight)
hy5edges[hy5edges$to %in% common,]

hist(data$edges$weight)
