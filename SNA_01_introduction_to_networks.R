## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 
## Social network analysis from positioned telemetry data
## COST action ETN Training School - České Budějovice (Česká republika)
## 
## PART 1 - INTRODUCTION TO NETWORKS IN R ####
##
## Author: Eneko Aspillaga (IMEDEA, CSIC-UIB)
## Contact: aspillaga@imedea.uib-csic.es
## Date: 18-20 May 2022
## 
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Set working directory (directory of the current script)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install required libraries:
# install.packages(c("igraph", "data.table", "lubridate", "sp", "sf", 
#                    "abind", "randomcoloR", "ggplot2", "ggraph"))

# Load libraries
library(igraph)
library(data.table)

# Function to rescale variables
rescale <- function(values, min, max) {
  range <- c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
  min + (values - range[1]) * (max - min) / (range[2] - range[1]) 
}


# 1. Download example data set =================================================

# The dataset we will use for this first session corresponds to the
# interactions between the characters of the 1st season of "Game of Thrones" 
# HBO TV Series (Beveridge and Chemers, 2018).
# Data source, details on how interactions were calculated, and a complete 
# analysis of all the seasons and book networks can be found in:
#     - https://networkofthrones.wordpress.com/
#     - https://github.com/mathbeveridge/gameofthrones
# Citation: Andrew Beveridge and Michael Chemers, "The Game of 'The Game of 
#           Thrones': Networked Concordances and Fractal Dramaturgy", 
#           in: Paola Brembilla and Ilaria De Pacalis (eds.), Reading 
#           Contemporary Serial Television Universes: A Narrative Ecosystem 
#           Framework, Routledge, 2018.

# Direct download from github
# nodes_path <- "https://raw.githubusercontent.com/mathbeveridge/gameofthrones/master/data/got-s1-nodes.csv"
# edges_path <- "https://raw.githubusercontent.com/mathbeveridge/gameofthrones/master/data/got-s1-edges.csv"
# nodes <- fread(nodes_path)
# edges <- fread(edges_path)

# Export to load them offline
# dir.create("./data/got_season1/", recursive = TRUE, showWarnings = FALSE)
# fwrite(nodes, "./data/got_season1/got-s1-nodes.csv")
# fwrite(edges, "./data/got_season1/got-s1-edges.csv")

# Load them from the exported files
nodes <- fread("./data/got_season1/got-s1-nodes.csv")
edges <- fread("./data/got_season1/got-s1-edges.csv")


# The "edges" data frame contains the number of interactions ("Weight") between 
# pairs of characters ("Source" and "Weight") in Season 1 of GoT. Interactions
# are undirected, so "Source" and "Target" variables are interchangeable.
head(edges) 

# The "nodes" data frame contains all the names of the characters in the
# network. Each character is represented by a code ("Id"), which is used
# in the "edges" dataset, and its name ("Label).
head(nodes)


# 2. Creating our first network with igraph ====================================

# Since our data is in a "data.frame" format, we will use the
# "graph_from_data_frame" function. The first two columns should correspond to 
# the nodes. The rest of variables are included as edge attributes. However, we 
# will change the spelling of the "Weight" column to "weight", so igraph 
# automatically understands it as edge weights.
colnames(edges)[colnames(edges) == "Weight"] <- "weight"
colnames(edges)
net <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
class(net)

# Extract nodes/vertices and check their attributes
V(net) 
vertex_attr_names(net)
vertex_attr(net)

# Extract links/edges and check their attributes
E(net)
edge_attr_names(net)
edge_attr(net)

# With the "V" and "E" functions, we can access directly to the attributes or 
# create new ones (important for plots)
V(net)$name
V(net)$Label
V(net)$extra_attribute <- 1:length(V(net))
vertex_attr_names(net)
V(net)$extra_attribute
E(net)$weight

# We can get back a data.frame with the edges
df <- as_data_frame(net)
head(df)

# Or an adjacency matrix of the network
adj <- as_adjacency_matrix(net, attr = "weight")
adj <- as.matrix(adj)
adj[1:10, 1:10]
isSymmetric(adj)


# 3. Extracting network- and node-level parameters =============================

# Network-level parameters: Useful to compare different networks

# Total number of edges:
length(E(net))
ecount(net)

# Edge-weight distribution:
quantile(E(net)$weight)
median(E(net)$weight)
hist(E(net)$weight)

# Edge-density:
edge_density(net)


# Node level parameters:

# Degree: number of adjacent edges (number of links of each node)
deg <- degree(net)
deg
hist(deg)

# Strength (weighted degree): sum of edge weights (number of interactions of
# each node)
str <- strength(net)
str
hist(str)

# Betweeness centrality: number of shortest paths going through the node.
# Indicates the importance of each node in connecting the network
bet <- betweenness(net)
bet
hist(bet)

# Eigenvector centrality: Number of important nodes to which each node is 
# connected (centrality of each node is proportional to the sum of the 
# centralities of the adjacent nodes)
eig <- eigen_centrality(net)
eig
eig$vector
hist(eig$vector)

# Page rank: Google Page Rank algorithm. Returns the importance of each node
# a percentage. The value for each node is proportional to the importance of
# all the nodes connected to it
ran <- page_rank(net)
ran
ran$vector
sum(ran$vector)
hist(ran$vector)


# Plot the characters with the highest values of each parameter
par(mar = c(3.1, 8.1, 3.1, 1), mfrow = c(2, 3))
barplot(deg[order(deg, decreasing = TRUE)][1:5], las = 2, horiz = TRUE, 
        main = "Degree")
barplot(str[order(str, decreasing = TRUE)][1:5], las = 2, horiz = TRUE,
        main = "Weighted degree")
barplot(bet[order(bet, decreasing = TRUE)][1:5], las = 2, horiz = TRUE,
        main = "Betweenness")
barplot(eig$vector[order(eig$vector, decreasing = TRUE)][1:5], las = 2, 
        horiz = TRUE, main = "Eigenvector centrality")
barplot(ran$vector[order(ran$vector, decreasing = TRUE)][1:5], las = 2, 
        horiz = TRUE, main = "Page Rank")


# 4. Plotting networks using the "igraph" package ==============================

# Using the "igraph" package:
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
plot(net)

# Plotting parameters can be modified by providing specific arguments to the 
# "plot" function or by adding new names attributes to the "igraph" object.
# Note that the arguments to modify vertices or edges are preceded by "vertex."
# or "edge." prefixes, respectively. More details in "help(igraph.plotting)".

# Modify node size by providing "vertex.size" argument. We will make it
# proportional to the betweenness of each node. We will also remove the vertex
# labels.
plot(net, vertex.size = rescale(log(bet+1), 0.1, 6), vertex.label = "")

# Same as before but adding the "size" attribute to vertices:
V(net)$size <- rescale(log(bet+1), 0.1, 6)
plot(net, vertex.label = "")

# Add vertex labels. We will use the name of each character
V(net)$label <- nodes$Label[match(V(net)$name, nodes$Id)]
plot(net)

# Modify the size of each label (label.cex) depending on their strength.
# When plotting, we will also change their font family and the color.
V(net)$label.cex <- rescale(log(str), 0.3, 1.8)
plot(net, vertex.label.family = "sans", vertex.label.color = "gray30")

# Modify edge with to make it proportional to its weight.
E(net)$width <- rescale(log(E(net)$weight), 0.5, 8)
plot(net, vertex.label = "")


# Change layout. There are different algorithms to place the nodes. They can
# be checked in "help(layout_)". Edges can be curved using the "edge.curved"
# argument.
plot(net, layout = layout_nicely, 
     edge.curved = 0.25, # Adds 25% of curvature to edges
     edge.color = adjustcolor("gray20", 0.6),
     vertex.label.family = "sans", vertex.label.color = "gray20",
     vertex.color = "coral")


# Searching for communities. We can look for internal structures within the
# network (groups of nodes related to each other) by applying different 
# algorithms. More information at help(communities).
groups <- cluster_fast_greedy(net)
class(groups)
plot_dendrogram(groups)

# We can extract the group to which belongs each node
membership(groups)
groups$membership
table(groups$membership)

# Add groups to the plot
plot(groups, net, layout = layout_nicely, 
     edge.curved = 0.25, # Adds 25% of curvature to edges
     edge.color = adjustcolor("gray20", 0.6),
     vertex.label = "")


# Change vertex color depending on the group
length(unique(groups$membership))
gr_col <- c("gray50", "coral", "gold", "yellowgreen", "dodgerblue", 
            "forestgreen", "darkorchid")
V(net)$color <- gr_col[groups$membership]

plot(net, layout = layout_with_graphopt, 
     edge.curved = 0.25, edge.color = adjustcolor("gray20", 0.6),
     vertex.label = "")


# Export the plot
jpeg("got_s1_network_igraph.jpeg", width = 4000, height = 4000, 
     res = 200)
set.seed(27)
plot(net, layout = layout_with_graphopt, 
     edge.curved = 0.25, 
     edge.color = adjustcolor("gray50", 0.6),
     vertex.label.family = "sans", vertex.label.color = "black",
     vertex.frame.color = "gray30")
dev.off()


# 5. Plotting networks using the "ggraph" package ==============================

# "ggraph" is an extension of the "ggplot2" package.
library(ggplot2)
library(ggraph)

# An nice guide to network plots with "ggraph" can be found at:
# https://mr.schochastics.net/material/netvizr/

# Basic plot
ggraph(net, layout = "stress") +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()

# Change plot parameters
ggraph(net, layout = "stress") + 
  geom_edge_link(aes(edge_width = weight), edge_colour = "gray50", ) +
  geom_node_point(aes(size = deg, fill = factor(groups$membership)), 
                  shape = 21) +
  geom_node_text(aes(label = label, size = deg), colour = "black", 
                 repel = TRUE) +
  scale_edge_width_continuous(range = c(0.2, 2)) +
  scale_size(range = c(0.5, 5))+
  scale_fill_brewer(palette = "Set2") +
  theme_graph() +
  theme(legend.position = "none")


# Using centrality layout. A centrality index is given ("cent") to create a 
# concentric layout where the most central nodes are put in the center
ggraph(net, layout = "centrality",  cent = eig$vector) + 
  geom_edge_link(aes(edge_width = weight), edge_colour = "gray50") +
  geom_node_point(aes(size = deg, fill = factor(groups$membership)), 
                  shape = 21) +
  geom_node_text(aes(label = label, size = deg), colour = "black", 
                 repel = TRUE) +
  scale_edge_width_continuous(range = c(0.2, 2)) +
  scale_size(range = c(0.5, 6))+
  scale_fill_brewer(palette = "Set2") +
  theme_graph() +
  theme(legend.position = "none")

