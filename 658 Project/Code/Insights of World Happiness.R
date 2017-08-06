# packages
install.packages("OpenMx")
install.packages("igraph")   # network constrcution
install.packages("RColorBrewer")  # nodes color

library(OpenMx)
library(scales)
library(igraph)
library(RColorBrewer)

mydata <- read.csv(file = "data.csv", header = TRUE)
cnames <- mydata$Country

# generate 6 matrix of diff of each feature
vec_to_matrix <- function(vector,names){
  veclen <- length(vector)
  mymat <- matrix(nrow = veclen, ncol = veclen)
  for (i in 1:veclen){
    for (j in 1:veclen){
        mymat[j,i] <- abs(vector[i]-vector[j])
    }
  }
  mymat <- apply(mymat, MARGIN = 2, FUN = function(X) 10 - ceiling((X - min(X))*10/diff(range(X))))
  row.names(mymat) <- names
  colnames(mymat) <- names
  return(mymat)
}

vec_to_matrix2 <- function(vector,names){
  veclen <- length(vector)
  mymat <- matrix(nrow = veclen, ncol = veclen)
  for (i in 1:veclen){
    for (j in 1:veclen){
      mymat[j,i] <- abs(vector[i]-vector[j])
    }
    mymat[,i] <- 1/(scale(mymat[,i])+3)
  }
  row.names(mymat) <- names
  colnames(mymat) <- names
  return(mymat)
}


# Binary edge detection for each matrix

Nth_large <- function(x, N){
  N <- min(N, length(x))
  x <- sort(x[x >= -sort(-x, partial=N)[N]][1:N])
  x[1]
}

edge_detect <- function(matrix,N){
  for (i in 1:ncol(matrix)){
    for (j in 1:nrow(matrix)){
      if (matrix[j,i] <= Nth_large(matrix[,i],N)){
        matrix[j,i] <- 0
      }
    }
  }
  return(matrix)
}

edge_detect2 <- function(matrix,N){
  for (i in 1:ncol(matrix)){
    for (j in 1:nrow(matrix)){
      if (matrix[j,i] <= Nth_large(matrix[,i],N)){
        matrix[j,i] <- 0
      }
      else{
        matrix[j,i] <- 1
      }
    }
  }
  return(matrix)
}
write.csv(adjmatrix,"adj2.csv")
matrix_economy <- vec_to_matrix(mydata$Economy..GDP.per.Capita.,cnames)
matrix_family <- vec_to_matrix(mydata$Family,cnames)
matrix_health <- vec_to_matrix(mydata$Health..Life.Expectancy.,cnames)
matrix_freedom <- vec_to_matrix(mydata$Freedom,cnames)
matrix_trust <- vec_to_matrix(mydata$Trust..Government.Corruption.,cnames)
matrix_generosity <- vec_to_matrix(mydata$Generosity,cnames)


adjmatrix <- (matrix_economy + matrix_family + matrix_health + matrix_freedom + matrix_trust + matrix_generosity)/10 

adjmatrix <- edge_detect(adjmatrix,6)

matrix_list <- c("matrix_economy","matrix_family","matrix_health","matrix_freedom","matrix_trust","matrix_generosity")

# set edges between nodes
for (i in (1:length(matrix_list))){
  assign(matrix_list[i],edge_detect2(get(matrix_list[i]),7))
}


# construct network of countries
net <- graph_from_adjacency_matrix(adjmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
net2 <- graph_from_adjacency_matrix(bigmatrix, mode = "undirected", weighted = TRUE, diag = FALSE)
net_freedom <- graph_from_adjacency_matrix(matrix_freedom, mode = "undirected", weighted = TRUE, diag = FALSE)

# community detection of net
greedy_comm = fastgreedy.community(net)
plot(as.dendrogram(greedy_comm))
V(net)$greedy_comm <- cut_at(greedy_comm,4)

between_comm = edge.betweenness.community(net, directed=F)
plot(as.dendrogram(between_comm),cex.color="blue")
cut_at(between_comm,4)
V(net)$between_comm <- cut_at(between_comm,4)

K_cluster <- kmeans(mydata[7:13], centers=4, nstart=10,iter.max = 100)
V(net)$Kcluster <- K_cluster$cluster

# three type of clustering
mydata$greedy_comm <- V(net)$greedy_comm
mydata$between_comm <- V(net)$between_comm
mydata$K_cluster <- V(net)$Kcluster

write.csv(mydata, file = "group_new.csv")

#community detection of net_freedom
greedy_comm_eco = fastgreedy.community(net_freedom)
plot(as.dendrogram(greedy_comm_eco))
V(net_freedom)$greedy_comm <- cut_at(greedy_comm,4)

# set ohter parameters of the net
V(net)$region <- mydata$Region
V(net)$size <- log(mydata$Population)-10
V(net)$lat <- mydata$Latitute
V(net)$lng <- mydata$longtitute
V(net)$GPP <- mydata$GPP/10000+3

# set other parameters of net_eco
V(net_freedom)$region <- mydata$Region
V(net_freedom)$size <- log(mydata$Population)-10
V(net_freedom)$lat <- mydata$Latitute
V(net_freedom)$lng <- mydata$longtitute
V(net_freedom)$GPP <- mydata$GPP/10000+3

# network visualization
colrs <- brewer.pal(4, "Spectral")
col <- adjustcolor(colrs, alpha=.6)


comm.graph <-  contract.vertices(net, V(net)$region, vertex.attr.comb=list(size="sum",name="ignore"))
V(comm.graph)$name <- region
plot(comm.graph,vertex.size=5)
comm.graph = simplify(comm.graph)
plot(comm.graph)

# igraph with world map
lo <- layout.norm(as.matrix(mydata[,14:15]),xmin = min(mydata[,14]), xmax = max(mydata[,14]), ymin = min(mydata[,15]), ymax = max(mydata[,15]))

plot(net, layout=lo,
     xlim = c(-180, 180), 
     ylim = c(-90, 90), 
     rescale = FALSE,
     vertex.size=V(net)$GPP, vertex.color=color[model1$cluster],
     vertex.label.font=1, vertex.label.cex=.5,
     edge.color="grey80")

# heatmap
netm <- as_adjacency_matrix(comm.graph, attr="weight",sparse = F)
colnames(netm) <- region
rownames(netm) <- region
palf <- colorRampPalette(c("gold","dark orange"))(5)
heatmap(as.matrix(netm), Rowv = NA, Colv = NA, col = palf, scale="none", margins=c(10,10) )

region <- c("Eastern Asia","North America","Southeastern Asia","Western Europe","Australia and New Zealand","Latin America and Caribbean","Sub-Saharan Africa","Middle East and Northern Africa","Central and Eastern Europe","Southern Asia")
region <- rev(region)

saveAsGEXF(net,"world.gexf")
saveAsGEXF(comm.graph,"my_graph_comm.gexf")
saveAsGEXF(net_freedom,"net_free.gexf")

aggregate(as.numeric(mydata$Population), by=list(mydata$between_comm), FUN=sum)
aggregate(mydata[8:13], by=list(mydata$between_comm), FUN=mean)
View(aggregate(mydata[8:13], by=list(mydata$Region), FUN=mean))
aggregate(as.numeric(mydata$Population),by=list(mydata$Region), FUN=sum)

# Calculate the mean of each cluster
grouped <- read.csv(file = "group.csv", header = TRUE)

# using fast-greedy result
features = c(grouped$Economy..GDP.per.Capita., grouped$Family, grouped$Health..Life.Expectancy., 
             grouped$Freedom, grouped$Trust..Government.Corruption., grouped$Generosity)

group1 <- subset(grouped, greedy_comm == 1, select= Economy..GDP.per.Capita. : Generosity)
group2 <- subset(grouped, greedy_comm == 2, select= Economy..GDP.per.Capita. : Generosity)
group3 <- subset(grouped, greedy_comm == 3, select= Economy..GDP.per.Capita. : Generosity)
group4 <- subset(grouped, greedy_comm == 4, select= Economy..GDP.per.Capita. : Generosity)

avg2 <- colMeans(group2, na.rm = FALSE, dims = 1)

# radar chart
summary(grouped)
maxmin <- data.frame(
  Economy=c(1.5, 0),  Family=c(1.5, 0),  Health=c(1.0, 0),
  Freedom=c(0.6, 0),  Trust=c(0.6, 0),  Generosity=c(0.6, 0))

rad1 = rbind(maxmin,colMeans(group1, na.rm = FALSE, dims = 1))
rad2 = rbind(maxmin,colMeans(group2, na.rm = FALSE, dims = 1))
rad3 = rbind(maxmin,colMeans(group3, na.rm = FALSE, dims = 1))
rad4 = rbind(maxmin,colMeans(group4, na.rm = FALSE, dims = 1))

install.packages('fmsb')
library(fmsb)
radarchart(rad1, axistype=0, seg=6, centerzero = TRUE, title = 'Group1')

# Custom the radarChart
par(mfrow=c(2, 2))
# Group 1 purple
radarchart( rad1 , axistype=1 , title = 'Group1',
            #custom polygon
            pcol=rgb(0.82,0.48,0.8,0.9) , pfcol=rgb(0.82,0.48,0.8,0.6) , plwd=4 , 
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.3,1.5,0.3), cglwd=0.8,
            #custom labels
            vlcex=0.8 )
# Group 2 orange
radarchart( rad2 , axistype=1, title = 'Group2', vlcex=0.8,
            pcol=rgb(1,0.45,0.24,0.9) , pfcol=rgb(1,0.45,0.24,0.6) , plwd=4 , 
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.3,1.5,0.3), cglwd=0.8)
# Group 3 green
radarchart( rad3 , axistype=1 , title = 'Group3', vlcex=0.8,
            pcol=rgb(0.41,0.68,0.21, 0.9) , pfcol=rgb(0.41,0.68,0.21,0.6) , plwd=4 , 
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.3,1.5,0.3), cglwd=0.8)
# Group 4 blue
radarchart( rad4 , axistype=1 , title = 'Group4', vlcex=0.8,
            pcol=rgb(0,0.7,0.85,0.9) , pfcol=rgb(0,0.7,0.85,0.6) , plwd=4 , 
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.3,1.5,0.3), cglwd=0.8)

dev.off()
