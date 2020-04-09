################ essais dap seq et connection a la base de donnees
load("D:/These/NetworkShiny/NetworkData/CO2DEGenes_faibleNitrate_CO2-N.RData")
setwd("D:/These/NetworkInference")

source("Funtions/Network_functions.R")

plotNetwork(data)

library(DBI)
library(RSQLite)
library(dplyr)
library(dbplyr)
library(ggplot2)

con <- dbConnect(RSQLite::SQLite(), "D:/These/DAPSeqData/databaseDapSeq.db")
src_dbi(con)


tfs <- data.frame(tbl(con, sql("SELECT * FROM tfs")))



query <- dbSendQuery(con, "SELECT COUNT(*), tf_agi, genes.name AS name, families.name AS familyName
                     FROM binding 
                     INNER JOIN genes ON binding.tf_agi=genes.agi 
                     INNER JOIN tfs ON binding.tf_agi=tfs.agi
                     INNER JOIN families ON tfs.family_id=families.id
                     WHERE amplified = 0 
                     GROUP BY tf_agi")
dbBind(query)
res <- data.frame(dbFetch(query))
res
res$x <- "x"

ggplotly(ggplot(res, aes(y=COUNT..., x=x, fill=x) )+ geom_violin( alpha=0.4)+ 
           geom_jitter(alpha=0.5, cex= 10,aes(fill=familyName))+
           theme(legend.position='none'), tooltip = c("name", "familyName"))




DapSeqApproved <- function(tf, target){
  query <- dbSendQuery(con, "SELECT * FROM binding WHERE tf_agi = ? AND gene_agi = ?")
  dbBind(query, list(tf, target))
  res <- data.frame(dbFetch(query))
  #print(res)
  if(dim(res)[1] > 0) return(TRUE)
  else(return(FALSE))
}


tf <- data$edges$from[1]
target <- data$edges$to[1]

DapSeqApproved(tf, target)



testedByDapSeq <- function(tf){
  query <- dbSendQuery(con, "SELECT * FROM tfs WHERE tfs.agi = ?")
  dbBind(query, list(tf))
  res <- data.frame(dbFetch(query))
  if(dim(res)[1] == 0) return(FALSE)
  else(return(TRUE))
}


data$edges$testedByDapSeq <- sapply(data$edges$from, testedByDapSeq)
sum(data$edges$testedByDapSeq)/dim(data$edges)[1]



data$edges$DapSeqAproved = mapply(DapSeqApproved, data$edges$from, data$edges$to)
sum(data$edges$DapSeqAproved) / sum(data$edges$testedByDapSeq)


dataRand <- data
dataRand$edges$to <- sample(dataRand$edges$to, replace=F)
dataRand$edges$DapSeqAproved = mapply(DapSeqApproved, dataRand$edges$from, dataRand$edges$to)



dataRandTested <- dataRand
dataRandTested$edges <- dataRandTested$edges[dataRandTested$edges$testedByDapSeq==T,]


ratesRandom = c()
for(i in 1:20){
  dataRandTested$edges <- dataRandTested$edges[dataRandTested$edges$testedByDapSeq==T,]
  dataRandTested$edges$to <- sample(dataRandTested$edges$to, replace=F)
  dataRandTested$edges$DapSeqAproved = mapply(DapSeqApproved, dataRandTested$edges$from, dataRandTested$edges$to)
  ratesRandom <- c(ratesRandom, sum(dataRandTested$edges$DapSeqAproved)/dim(dataRandTested$edges)[1])
}

df <- data.frame(rates = ratesRandom)

ggplot(df, aes(x=ratesRandom)) + geom_density(fill = "red", alpha=0.2)+ 
  geom_vline(xintercept  = sum(data$edges$DapSeqAproved) / sum(data$edges$testedByDapSeq)) +
  ggtitle("DAP Seq approuved links placed on DAP Seq aproved random links (CO2*N)")


hist(ratesRandom)



sum(dataRand$edges$DapSeqAproved) / sum(dataRand$edges$testedByDapSeq)

#tested <- data$edges[data$edges$testedByDapSeq==TRUE,]
#validated <- mapply(DapSeqApproved, tested$from, tested$to)


colorEdge <- function(tested, validated){
  if(tested & validated) return("red")
  if(tested & !validated) return("black")
  if(!tested & !validated) return("white")
}

data$edges$color <- mapply(colorEdge,  data$edges$testedByDapSeq, data$edges$DapSeqAproved)
 #<- ifelse(data$edges$DapSeqAproved, "blue", "#EEEEEE")

data$nodes$label <- data$nodes$Ontology
visNetwork(nodes = data$nodes, edges = data$edges)%>% 
  #visEdges(smooth = FALSE, color = '#333366') %>% 
  visEdges(smooth = FALSE, arrows = 'to') %>% 
  visPhysics(solver = "forceAtlas2Based", timestep = 0.9, minVelocity=10, 
             maxVelocity = 10, stabilization = F)%>%
  visOptions(selectedBy = "group", highlightNearest = TRUE,nodesIdSelection  = TRUE, collapse = TRUE)%>% 
  visEvents(click = "function(nodes){
                  Shiny.onInputChange('click', nodes.nodes);
                  ;}") %>% 
  visGroups(groupname = "Regulator", size = 28,
            color = list("background" = "grey", "border"="#CCCCCC"), shape = "square") %>% 
  visGroups(groupname = "Target Gene", color = "lightgrey") %>% 
  visNodes(borderWidth=0.5, font=list("size"=36)) 

plotNetwork(data)


save(data, file = "D:/These/NetworkShiny/NetworkData/CO2DEGenes_IronStarv_CO2-Fe.RData")

data$edges[1:10,]
sum(validated)


# verifie que binding est bien renseignee
query <- dbSendQuery(con, "SELECT * FROM binding WHERE tf_agi = ? AND gene_agi = ?")
dbBind(query, list("AT3G56400", "AT1G01010"))
dbFetch(query)



# requetes stylees
query <- dbSendQuery(con, "SELECT * FROM tfs INNER JOIN families ON families.id=tfs.family_id INNER JOIN genes ON tfs.agi=genes.agi  WHERE tfs.agi = ?")
query <- dbSendQuery(con, "SELECT * FROM tfs INNER JOIN families ON families.id=tfs.family_id WHERE tfs.agi = ?")
dbBind(query, list(tf))
dbFetch(query)

dbDisconnect(con)


dim(data$edges)

# combien de cibles dap seq on a avec HY5
query <- dbSendQuery(con, "SELECT * FROM binding WHERE tf_agi = ? AND amplified = 1")
dbBind(query, list("AT5G11260"))

res <- data.frame(dbFetch(query))

length(intersect(res$gene_agi, hy5))
