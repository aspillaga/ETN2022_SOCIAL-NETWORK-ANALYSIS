## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 
## Social network analysis from positioned telemetry data
## COST action ETN Training School - České Budějovice (Česká republika)
## 
## PART 2 - CREATING NETWORKS FROM TRAJECTORY DATA ####
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
library(lubridate)
library(sp)
library(sf)
library(abind)
library(randomcoloR)

# Function to rescale variables
rescale <- function(values, min, max) {
  range <- c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
  min + (values - range[1]) * (max - min) / (range[2] - range[1]) 
}


# 1. Load trayectory data ======================================================

# The dataset corresponds to high-resolution trajectories of razorfish 
# (Xyrichtys novacula) individual tracked in the Palma Bay Marine Reserve
# (Mallorca, Balearic Islands, Spain). High-resolution trajectories were
# obtained using the JSATS acoustic telemetry system (Lotek Wireless Inc.). 
# Positions were estimated using the UMAP software (also from Lotek) and post-
# processed by applying a positioning filter and a continuous-time correlated
# random walk movement model. The example data set consists on 153  individual 
# trajectories of one-day duration and with an interval between positions of 1 
# minute. 

# Details on the high-resolution telemetry system and the used post-processing
# method can be found at Aspillaga et al. 2021a:
# https://doi.org/10.1186/s40317-020-00224-w

# Details on the data set used in this script and the SNA applied to the data
# can be found at Aspillaga et al. 2021b:
# https://doi.org/10.3389/fmars.2021.688010

# Load fish reference file
fish_ref <- fread("./data/jsats_tracks/jsats_fish_ref_20190517.csv")
fish_ref

# Load trajectory data
tracks <- fread("./data/jsats_tracks/jsats_tracks_20190517.csv")
tracks


# 2. Explore trajectories ======================================================

# There are 153 individuals in the data set, 92 females and 61 males
nrow(fish_ref)
table(fish_ref$sex)

# Check number of positions per fish. Individuals have between 665 and 895 
# positions
n_pos <- table(tracks$tag.id)
range(n_pos)
hist(n_pos)

# Check time interval between positions (1 min)
time_int <- sapply(fish_ref$tag.id, function(id) {
  diff <- diff(tracks$date.time[tracks$tag.id == id])
  return(range(diff))
})
time_int
range(time_int)

# Check temporal coverage
# Calculate nummber of individuals with positions in each 1-min interval
n_ind <- tapply(tracks$tag.id, tracks$date.time, function(x) length(unique(x)))
plot(ymd_hms(names(n_ind)), n_ind, type = "l", xlab = "Time", ylab = "N ind.")
abline(h = 130, col = "red")

# Select time intervals with data for more than 130 individuals (~95% of total
# individuals)
time_range <- range(ymd_hms(names(n_ind))[n_ind > 130])
time_range
tracks <- tracks[tracks$date.time >= time_range[1] & 
                   tracks$date.time <= time_range[2], ]


# Plot the trajectories
# Select colors for individuals
col_females <- randomColor(sum(fish_ref$sex == "F"), hue = "orange")
col_males <- randomColor(sum(fish_ref$sex == "M"), hue = "blue")
fish_ref$color <- ""
fish_ref$color[fish_ref$sex == "F"] <- col_females
fish_ref$color[fish_ref$sex == "M"] <- col_males

# Plot all the tracks
plot(y ~ x, data = tracks, type = "n", asp = 1, ann = FALSE)
for (i in 1:nrow(fish_ref)) {
  lines(y ~ x, data = tracks[tracks$tag.id == fish_ref$tag.id[i]], 
        col = fish_ref$color[i])
}

# Zoom
plot(y ~ x, data = tracks, type = "n", asp = 1, ann = FALSE,
     xlim = c(476110, 476190), ylim = c(4368620, 4368710))
for (i in 1:nrow(fish_ref)) {
  lines(y ~ x, data = tracks[tracks$tag.id == fish_ref$tag.id[i]], 
        col = fish_ref$color[i])
}


# 3. Calculate associations ====================================================

# We will estimate associations between individuals based on their proximity.
# An association will be assumed when the distance between two individuals is
# less than 3 m.

# Convert tracks to SpatialPointsDataFrame (easier to compute distances)
proj <- CRS("+proj=utm +datum=WGS84 +zone=31") # Coordinarte system
tracks <- SpatialPointsDataFrame(coords = tracks[, c("x", "y")], data = tracks, 
                                 proj4string = proj)

# Extract vector of times (time-steps for computing associations)
times <- sort(unique(tracks$date.time))

# Estimate the distance between all the pairs of individuals in each time-step
dist_list <- lapply(times, function(t) {
  
  # Show the time interval that is being processed
  cat(format(t, "%H:%M:%S"), "\r")
  
  # Subset data for each time interval
  tr_sub <- tracks[tracks$date.time == t, ]
  
  # Calculate distances between individuals
  dist_sub <- st_distance(st_as_sf(tr_sub))
  
  colnames(dist_sub) <- tr_sub$tag.id
  rownames(dist_sub) <- tr_sub$tag.id
  
  # Match the rows and the columns of the matrix to include all the individuals
  # in the data set (individuals in the "fish_ref" data frame)
  dist_sub2 <- dist_sub[match(fish_ref$tag.id, rownames(dist_sub)),
                        match(fish_ref$tag.id, colnames(dist_sub))]
  colnames(dist_sub2) <- fish_ref$tag.id
  rownames(dist_sub2) <- fish_ref$tag.id
  
  # Return matrix with distances
  return(dist_sub2)
})

# Number of matrices in the list
length(dist_list) # 824 1-min intervals
dim(dist_list[[1]]) # Matrices of 153x153 individuals
dist_list[[1]][1:10, 1:10] # NAs where introduced when one any of the
                           # individuals was not present in the dataset

# Merge all the matrices in a 3-D array
distances <- abind(dist_list, along = 3)
dim(distances)

# Association threshold. Association will be considered if the distance between
# two individuals is smaller than 3 m. Similar to the positioning error.
assoc_th <- 3
assoc <- copy(distances)
assoc[which(distances > assoc_th)] <- 0
assoc[which(distances <= assoc_th)] <- 1
assoc

# Calculate the number of associations and observations of each dyad.
# Associations are calculated as the number of 1-min intervals in which two
# individuals (a dyad) are at less than 3m from each other. Observations are
# calculated as the total number of 1-min intervals in which two individuals
# where observed within the acoustic array.
associations <- apply(assoc, c(1, 2), sum, na.rm = TRUE)
observations <- apply(assoc, c(1, 2), function(x) sum(!is.na(x)))

dim(associations)
dim(observations)

associations[1:10, 1:10]
observations[1:10, 1:10]

# Calculate the association indices by dividing the number of associations by
# the total number of observations
assoc_index <- associations / observations
range(assoc_index)
assoc_index[1:10, 1:10]

# Plot the adjacency matrix as a heatmap
image(log(assoc_index))


# 4. Create the network and estimate parameters ================================

# In this case, since our data is in an adjacency matrix, we will generate
# the "igraph" object using the "graph_from_adjacency_matrix" function
net <- graph_from_adjacency_matrix(assoc_index, mode = "undirected", 
                                   diag = FALSE, weighted = TRUE)

# Check vertex attributes
vertex_attr_names(net) # The only attribute is the fish ID ("name")
V(net)$name

# Check if the "name" attribute matches the IDs in the reference file 
# ("fish_ref")
V(net)$name == fish_ref$tag.id

# Add vertex attributes (two ways: directly or using the "vertex_attr" function)
V(net)$sex <- fish_ref$sex
vertex_attr(net, "length") <- fish_ref$length.cm

# Also, we will add a new edge attribute that indicates the type of edge
# depending on the sex of individuals ("F-F": association between two females;
# "M-M": association between two males; and "M-F": association between 
# individuals of different sex).

# First, the edges in data frame format
edges <- as_data_frame(net)
head(edges)
from_sex <- fish_ref$sex[match(edges$from, fish_ref$tag.id)]
to_sex <- fish_ref$sex[match(edges$to, fish_ref$tag.id)]
E(net)$type <- ""
E(net)$type[from_sex != to_sex] <- "M-F"
E(net)$type[from_sex == "F" & to_sex == "F"] <- "F-F"
E(net)$type[from_sex == "M" & to_sex == "M"] <- "M-M"
E(net)$type

# Extract degree and strength of associations
deg <- degree(net)
str <- strength(net)

# Plots of degree and strength depending on the sex
boxplot(deg ~ V(net)$sex, outline = FALSE)
boxplot(str ~ V(net)$sex, outline = FALSE)

# Network parameters. First extract all the edges
edges <- as_data_frame(net)

# Number of edges depending on the type
barplot(table(edges$type), ylab = "No. of edges")

# Weigth of the edges depending on the type
boxplot(log(edges$weight) ~ edges$type)

# Separate networks depending on the type
net_mm <- delete_vertices(net, V(net)$sex == "F")
V(net_mm)$sex
unique(E(net_mm)$type)
net_ff <- delete_vertices(net, V(net)$sex == "M")
V(net_ff)$sex
unique(E(net_ff)$type)
net_mf <- delete_edges(net, E(net)[E(net)$type %in% c("M-M", "F-F")])
degree(net_mf)
unique(E(net_mf)$type)

# Calculate node parameters for each sub-network
deg_mm <- degree(net_mm)
deg_ff <- degree(net_ff)
deg_mf <- degree(net_mf)[V(net_mf)$sex == "M"]
deg_fm <- degree(net_mf)[V(net_mf)$sex == "F"]
boxplot(deg_mm, deg_ff, deg_mf, deg_fm, names = c("M-M", "F-F", "M-F", "F-M"),
        ylab = "Degree", outline = FALSE)

str_mm <- strength(net_mm)
str_ff <- strength(net_ff)
str_mf <- strength(net_mf)[V(net_mf)$sex == "M"]
str_fm <- strength(net_mf)[V(net_mf)$sex == "F"]
boxplot(str_mm, str_ff, str_mf, str_fm, names = c("M-M", "F-F", "M-F", "F-M"),
        ylab = "Strength", outline = FALSE)


# 5. Plot the network ==========================================================

# Node size as function of the size of the fish
V(net)$size <- rescale(V(net)$length, 0.4, 10)
plot(net, layout = layout_nicely, vertex.label = "", edge.curved = 0.25,
     edge.color = adjustcolor("gray30", 0.5))

# Add color for each fish
V(net)$color <- fish_ref$color
plot(net, layout = layout_nicely, vertex.label = "", edge.curved = 0.25,
     edge.color = adjustcolor("gray30", 0.5))

# Modify node shape depending on the sex
V(net)$shape <- ifelse(V(net)$sex == "F", "circle", "square")
plot(net, layout = layout_nicely, vertex.label = "", edge.curved = 0.25,
     edge.color = adjustcolor("gray30", 0.5))

# Width of the edge as function of its weight
E(net)$width <- rescale(log(E(net)$weight), 0.1, 4)
plot(net, layout = layout_nicely, vertex.label = "", edge.curved = 0.25,
     edge.color = adjustcolor("gray30", 0.5))

# Group individuals
groups <- cluster_fast_greedy(net, weights = E(net)$weight)
plot(groups, net, layout = layout_nicely, vertex.label = "", 
     edge.curved = 0.25)


# Simplify the network to show the strongest associations

# Remove individuals with no interactions
net_sub <-  delete_vertices(net, degree(net) == 0)
groups_sub <- cluster_fast_greedy(net_sub)

plot(groups_sub, net_sub, layout = layout_nicely, vertex.label = "", 
     edge.curved = 0.25)

# Show the 25% of the strongest edges
w_th <- quantile(E(net)$weight, 0.75)
w_th
net_sub2 <- delete_edges(net, E(net)[E(net)$weight < w_th])
net_sub2 <-  delete_vertices(net_sub2, degree(net_sub2) == 0)
groups_sub2 <- cluster_fast_greedy(net_sub2)
E(net_sub2)$width <- rescale(log(E(net_sub2)$weight), 0.1, 4)

plot(groups_sub2, net_sub2, layout = layout_nicely, vertex.label = "", 
     edge.curved = 0.5)


# Export data for to be used in the next script ################################
save(tracks, fish_ref, net, file = "./data/script_02_output.rda")

