load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")
setwd("D:/These/NetworkInference")

library(visNetwork)
library(igraph)
source("Funtions/Network_functions.R")


data$edges$color = "grey"
plotNetwork(data)


see_backbone <- function(file){
  back <- read.table(file, sep = '\t', h = T)
  colnames(back)[1:2] <- c("from", "to")
  head(back)
  data$edges$color <- ifelse(paste(data$edges$from, data$edges$to) %in% paste(back$from, back$to), "darkred", "#CCCCCC")
  plotNetwork(data)
  
}

see_backbone("D:/These/Backboning/networkCO2_N_3000.txt_2.32_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N.txt1.28_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N_top0.05.txt1.28_nc.csv")
see_backbone("D:/These/Backboning/NetworkCO2_N_top0.05.txt2.32_nc.csv")

