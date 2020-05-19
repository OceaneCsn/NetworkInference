load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")
setwd("D:/These/NetworkInference")

library(visNetwork)
library(igraph)
source("Funtions/Network_functions.R")


data$edges$color = "grey"
plotNetwork(data)



see_backbone <- function(file) {
  back <- read.table(file, sep = '\t', header = T)
  colnames(back)[1:2] <- c("from", "to")
  head(back)
  data$edges$color <-
    ifelse(
      paste(data$edges$from, data$edges$to) %in% paste(back$from, back$to),
      "darkred",
      "#CCCCCC"
    )
  plotNetwork(data)
  
}



file = "D:/These/Backboning/networkCO2_N_3000.txt_2.32_nc.csv"
backboneData <- function(file) {
  back <- read.table(file, sep = '\t', header = T)
  colnames(back)[1:2] <- c("from", "to")
  newData <- data
  newData$edges <-
    newData$edges[paste(data$edges$from, data$edges$to) %in% paste(back$from, back$to), ]
  newData$nodes <-
    newData$nodes[newData$nodes$id %in% unique(union(newData$edges$from, newData$edges$to)), ]
  print('edges :')
  print(dim(newData$edges))
  print('nodes :')
  print(dim(newData$nodes))
  return(newData)
}



newData <- backboneData(file)
data <- newData
save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N_Backbone_p0.01.RData")



see_backbone("D:/These/Backboning/networkCO2_N_3000.txt_2.32_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N.txt1.28_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N_top0.05.txt1.28_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N_top0.05.txt2.32_nc.csv")


####################################### network with TFs only

newData <- data

newData$nodes <- newData$nodes[newData$nodes$group == "Regulator", ]
dim(newData$nodes)
newData$edges <-
  newData$edges[newData$edges$from %in% newData$nodes$id &
                  newData$edges$to %in% newData$nodes$id , ]

newData$edges$color = "grey"
plotNetwork(newData)

data <- newData
save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N_TFs_Only.RData")




  