### IGraph
http://kateto.net/networks-r-igraph
### With Igraph
library(igraph)

AMR_Conv_samples = which(pData(AMR_raw_analytic_data[[1]])$Treatment == "Conventional")
AMR_RWA_samples = which(pData(AMR_raw_analytic_data[[1]])$Treatment == "Natural")

AMR_class_Conv <- AMR_raw_analytic_data[[1]][, AMR_Conv_samples]
AMR_class_RWA <- AMR_raw_analytic_data[[1]][, AMR_RWA_samples]


Mech_drx_net2 <- graph.incidence(t(MRcounts(AMR_class_RWA)),weighted = TRUE) #this is a weighted network
class_tax <- fData(AMR_analytic_data[[1]])
class_tax$ID <- class_tax$Feature
class_tax$short_ID <- class_tax$Feature


nodes <- rbind(metadata, class_tax[,-1],fill=TRUE)

#net <- graph_from_data_frame(d=E(Mech_drx_net2), vertices=nodes, directed=T) 


plot(Mech_drx_net2, vertex.label.color="black", vertex.label.dist=1,vertex.size=7,layout=layout.auto)

# Set node size based on audience size:
V(Mech_drx_net2)$size <- V(Mech_drx_net2)$audience.size*0.7

## bipartite network
Mech_drx_net2.bp <- bipartite.projection(Mech_drx_net2)
print(Mech_drx_net2.bp[[1]], g=TRUE, e=TRUE)

plot(Mech_drx_net2.bp$proj1, vertex.label.color="black", vertex.label.dist=1,vertex.size=7)
plot(Mech_drx_net2.bp$proj2, vertex.label.color="black", vertex.label.dist=1,vertex.size=7,layout=layout_nicely)
plot(Mech_drx_net2.bp$proj2, vertex.label.color="black", vertex.label.dist=1,vertex.size=7,layout=layout.auto)

cut.off <- mean(Mech_drx_net2$weight) 

net.sp <- delete_edges(net, E(net)[weight<cut.off])

plot(Full.table) 

V(Mech_drx_net2)$name
V(Mech_drx_net2)$label.color
V(Mech_drx_net2)$type
E(Mech_drx_net2)$weight
E(Mech_drx_net2)$type

plot(Mech_drx_net2, edge.arrow.size=.4)
Mech_drx_net2 <- simplify(Mech_drx_net2, remove.multiple = F, remove.loops = T) 
plot(Mech_drx_net2, edge.arrow.size=.4,vertex.label=NA)


plot(Mech_drx_net2, layout=-layout.bipartite(Mech_drx_net2)[,2:1], vertex.size=30, vertex.shape=ifelse(V(Mech_drx_net2)$type, "rectangle", "circle"), vertex.color=ifelse(V(Mech_drx_net2)$type, "red", "cyan"))
plot(Mech_drx_net2, layout= layout.reingold.tilford(Mech_drx_net2,circular=TRUE), vertex.size=20, vertex.color="yellow", edge.width=E(Mech_drx_net2)$weight)

Mech_drx_net2_fr <- layout_with_fr(Mech_drx_net2)
plot(Mech_drx_net2, layout=Mech_drx_net2_fr )

hist(E(Mech_drx_net2)$weight)
mean(E(Mech_drx_net2)$weight)
sd(E(Mech_drx_net2)$weight)

cut.off <- mean(E(Mech_drx_net2)$weight)
net.sp <- delete_edges(Mech_drx_net2, E(Mech_drx_net2)[weight<cut.off])
plot(net.sp) 

hs <- hub_score(Mech_drx_net2, weights=NA)$vector
as <- authority_score(Mech_drx_net2, weights=NA)$vector
par(mfrow=c(1,2))
plot(Mech_drx_net2, vertex.size=hs*50, main="Hubs")
plot(Mech_drx_net2, vertex.size=as*30, main="Authorities")

### Community detection
ceb <- cluster_edge_betweenness(Mech_drx_net2) 
dendPlot(ceb, mode="hclust")
plot(ceb,Mech_drx_net2)
class(ceb)
membership(ceb)

##clustering based on labels
clp <- cluster_label_prop(Mech_drx_net2)
plot(clp, Mech_drx_net2)

## based on greedy optimization of modularity
cfg <- cluster_fast_greedy(as.undirected(Mech_drx_net2))
plot(cfg, as.undirected(Mech_drx_net2))

## K core decomposition
kc <- coreness(Mech_drx_net2, mode="all")
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

plot(Mech_drx_net2, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])

cluster <- fastgreedy.community(Mech_drx_net2)
length(cluster)
sizes(cluster)
membership(cluster)
plot(cluster,Mech_drx_net2)
library(ape)
dendPlot(cluster,mode="phylo")
dendPlot(cluster,mode="hclust")

clust1 <- (cluster,Mech_drx_net2)
#write.graph(g=Mech_drx_net2,file= 'Mechanism_drx_network.gml', format='gml')
l <- layout.kamada.kawai(Mech_drx_net2)
plot(Mech_drx_net2, layout=l)


setwd("C:\\Users\\enrique\\Dropbox\\Projects\\NCBA2\\ManuscriptAnalysis\\Oct2016\\figures\\")

install.packages("sna","/home/enrique/anaconda3/envs/compute/lib/R/library",repos='http://cran.rstudio.com/')
A <- get.adjacency(Mech_drx_net2, sparse=FALSE)
library(network)
g <- network::as.network.matrix(A)
library(sna)
pdf("network_Resistome.pdf")
sna::gplot.target(g, degree(g), main="Degree", circ.col="skyblue",usearrows = FALSE,vertex.col=ifelse(V(Mech_drx_net2)$type, "red", "cyan"),edge.col="darkgray")
sna::gplot.target(g, evcent(g)$vector, main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Mech_drx_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
sna::gplot.target(g, betweenness(g), main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Mech_drx_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
sna::gplot.target(g, closeness(g), main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Mech_drx_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
dev.off()




