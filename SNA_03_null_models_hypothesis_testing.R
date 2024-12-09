## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 
## Social network analysis from positioned telemetry data
## COST action ETN Training School - České Budějovice (Česká republika)
## 
## PART 3 - NULL MODELS FOR HYPOTHESIS TESTING ####
##
## Author: Eneko Aspillaga (IMEDEA, CSIC-UIB)
## Contact: aspillaga@imedea.uib-csic.es
## Date: 18-20 May 2022
## 
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Set working directory (directory of the current script)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install required libraries:
# install.packages(c("igraph", "data.table", "lubridate", "sp", "rgeos", 
#                    "abind", "randomcoloR", "ggplot2", "ggraph", 
#                    "adehabitatHR"))

# Load libraries
library(igraph)
library(data.table)
library(lubridate)
library(sp)
library(rgeos)
library(abind)
library(adehabitatHR)


# Load data from the previous script
load("./data/script_02_output.rda")

# Function to rescale variables
rescale <- function(values, min, max) {
  range <- c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
  min + (values - range[1]) * (max - min) / (range[2] - range[1]) 
}


# 1. Generate random trajectories ==============================================

# Here, we will use a pre-network randomization approach. We will shuffle the
# location of all the trajectories to break the spatial connections but
# keeping their general structure. To do this, we will calculate the mean X and
# Y coordinates of each individual trajectories and displace it to match a 
# different means. 

# Calculate mean coordinates for each individual
fish_ref$mean_x <- NA
fish_ref$mean_y <- NA
for (i in 1:nrow(fish_ref)) {
  fish_ref$mean_x[i] <- mean(tracks$x[tracks$tag.id == fish_ref$tag.id[i]])
  fish_ref$mean_y[i] <- mean(tracks$y[tracks$tag.id == fish_ref$tag.id[i]])
}

fish_ref

# Plot mean coordinates of each individual
plot(fish_ref$mean_x, fish_ref$mean_y, col = fish_ref$color, pch = 16, 
     cex = rescale(fish_ref$length.cm, 0.5, 2), asp = 1)

# Calculate minimum convex polygon of all the mean coordinates. Later, we will
# restrict the sampling of the new coordinates to this area
mcp <- mcp(SpatialPoints(coords = fish_ref[, c("mean_x", "mean_y")]), 
           percent = 95)
plot(mcp, add = TRUE, border = "red")


# Define a function to randomize tracks. It will randomly sample a new pair
# of mean coordinates for each individual inside the MCP and then center the
# entire track
randomTracks <- function() {
  
  # Sample new coordinates
  new_coord <- spsample(mcp, n = nrow(fish_ref), type = "random")
  new_coord <- coordinates(new_coord)
  
  # Plot new coordinates
  # plot(new_coord, col = fish_ref$color, pch = 16, 
  #      cex = rescale(fish_ref$length.cm, 0.5, 2), asp = 1)
  
  # Re-center all the trajectories to the new coordinates
  new_tracks <- data.frame(tracks)
  for (i in 1:nrow(fish_ref)) {
    x_shift <- fish_ref$mean_x[i] - new_coord[i, "x"]
    y_shift <- fish_ref$mean_y[i] - new_coord[i, "y"]
    indx <- new_tracks$tag.id == fish_ref$tag.id[i]
    new_tracks$x[indx] <- new_tracks$x[indx] - x_shift
    new_tracks$y[indx] <- new_tracks$y[indx] - y_shift
  }
  return(new_tracks)
}

# Generate one iteration and plot original vs random tracks
set.seed(29)
random_tr <- randomTracks()

par(mfrow = c(1, 2))
plot(y ~ x, data = tracks, type = "n", asp = 1, main = "Original tracks")
for (i in 1:nrow(fish_ref)) {
  lines(y ~ x, data = tracks[tracks$tag.id == fish_ref$tag.id[i], ], 
        col = fish_ref$color[i])
}
plot(y ~ x, data = random_tr, type = "n", asp = 1, main = "Radomized tracks")
for (i in 1:nrow(fish_ref)) {
  lines(y ~ x, data = random_tr[random_tr$tag.id == fish_ref$tag.id[i], ], 
        col = fish_ref$color[i])
}


# Function to calculate an adjacency matrix of associations from tracks (same as 
# in the previous script)
adjMat <- function(tr) {
  
  # Convert tracks to SpatialPointsDataFrame (easier to compute distances)
  proj <- CRS("+proj=utm +datum=WGS84 +zone=31")
  tr <- SpatialPointsDataFrame(coords = tr[, c("x", "y")], 
                               data = tr, proj4string = proj)
  
  # Extract vector of times (time-steps for computing associations)
  times <- sort(unique(tr$date.time))
  
  # Estimate the distance between all the pairs of individuals in each time-step
  dist_list <- lapply(times, function(t) {
    
    # Show the time interval that is being processed
    cat(format(t, "%H:%M:%S"), "\r")
    
    # Subset data for each time interval
    tr_sub <- tr[tr$date.time == t, ]
    
    # Calculate distances between individuals
    dist_sub <- gDistance(tr_sub, byid = TRUE)
    colnames(dist_sub) <- tr_sub$tag.id
    rownames(dist_sub) <- tr_sub$tag.id
    
    # Match the rows and the columns of the matrix to include all the 
    # individuals in the data set (individuals in the "fish_ref" data frame)
    dist_sub2 <- dist_sub[match(fish_ref$tag.id, rownames(dist_sub)),
                          match(fish_ref$tag.id, colnames(dist_sub))]
    colnames(dist_sub2) <- fish_ref$tag.id
    rownames(dist_sub2) <- fish_ref$tag.id
    
    # Return matrix with distances
    return(dist_sub2)
  })
  
  # Merge all the matrices in a 3-D array
  distances <- abind(dist_list, along = 3)
  
  # Association threshold. Association will be considered if the distance 
  # between two individuals is smaller than 3 m
  assoc_th <- 3
  assoc <- copy(distances)
  assoc[which(distances > assoc_th)] <- 0
  assoc[which(distances <= assoc_th)] <- 1
  
  # Calculate the number of associations and observations of each dyad.
  # Associations are calculated as the number of 1-min intervals in which two
  # individuals (a dyad) are at less than 3m from each other. Observations are
  # calculated as the total number of 1-min intervals in which two individuals
  # where observed within the acoustic array.
  associations <- apply(assoc, c(1, 2), sum, na.rm = TRUE)
  observations <- apply(assoc, c(1, 2), function(x) sum(!is.na(x)))
  
  # Calculate the association indices by dividing the number of associations by
  # the total number of observations
  assoc_index <- associations / observations
  return(assoc_index)
}

# Example:
null_assoc <- adjMat(random_tr)
null_assoc[1:10, 1:10]


# Function to generate a network from an adjacency matrix and add the node and
# edge attributes (same as in the previous script)
setNetwork <- function(assoc_index) {
  n <- graph_from_adjacency_matrix(assoc_index, mode = "undirected", 
                                   diag = FALSE, weighted = TRUE)
  # Add node attributes
  V(n)$sex <- fish_ref$sex
  V(n)$length <- fish_ref$length.cm
  
  # Add "type" attribute to edges
  edges <- as_data_frame(n)
  head(edges)
  from_sex <- V(net)$sex[match(edges$from, V(net)$name)]
  to_sex <- V(net)$sex[match(edges$to,  V(net)$name)]
  E(n)$type <- ""
  E(n)$type[from_sex != to_sex] <- "M-F"
  E(n)$type[from_sex == "F" & to_sex == "F"] <- "F-F"
  E(n)$type[from_sex == "M" & to_sex == "M"] <- "M-M"
  
  return(n)
}


# Example
null_net <- setNetwork(null_assoc)
plot(null_net, vertex.label = "")

# We can now estimates the network parameters
degree(null_net)


# We will use "lapply" to generate a large number of random networks. The
# following code chunk generates 3 random networks and stores them in a list
null_net_list <- lapply(1:3, function(i) {
  cat("\nGeneraring null network no.", i, "\n")
  random_tr <- randomTracks()
  null_assoc <- adjMat(random_tr)
  null_net <- setNetwork(null_assoc)
  return(null_net)
})

# We can now estimates the network parameters of each network separately
degree(null_net_list[[1]])
degree(null_net_list[[2]])
degree(null_net_list[[3]])
boxplot(lapply(null_net_list, degree))

# Or process all of them at once using the "lapply" function
mean_deg <- lapply(null_net_list, function(n) {
  mean(degree(n))
})
mean_deg # Returns a list with single values
unlist(mean_deg) # Convert the list to a vector

# Plot the degree distribution using the "lapply" function
boxplot(lapply(null_net_list, degree))

# Apply the functions to generate 500 null-networks. This part takes quite a
# long time to run, soyou can directly load the resulting list with null 
# networks running the code some lines below.
# set.seed(28)
# null_net_list <- lapply(1:500, function(i) {
#   cat("\nGeneraring null network no.", i, "\n")
#   random_tr <- randomTracks()
#   null_assoc <- adjMat(random_tr)
#   null_net <- setNetwork(null_assoc)
#   return(null_net)
# })
# saveRDS(null_net_list, file = "./data/random_networks_n500.rds")

# Load the 500 random networks
null_net_list <- readRDS("./data/random_networks_n500.rds")
length(null_net_list)


# 2. Extract and compare null and observed parameters ==========================

# We will compare network parameters: No of edges, mean edge weight
# We will compare the node parameters: mean degree, mean strength

# Function to extract all the parameters from the network:

getParam <- function(net) {
  
  # Separate networks depending on the type
  net_mm <- delete.vertices(net, V(net)$sex == "F")
  net_ff <- delete.vertices(net, V(net)$sex == "M")
  net_mf <- delete.edges(net, E(net)[E(net)$type %in% c("M-M", "F-F")])
  
  # Network parameters:
  # Number of edges
  n_edges_mm <- length(E(net_mm))
  n_edges_ff <- length(E(net_ff))
  n_edges_mf <- length(E(net_mf))
  
  # Mean weight
  w_edges_mm <- mean(E(net_mm)$weight)
  w_edges_ff <- mean(E(net_ff)$weight)
  w_edges_mf <- mean(E(net_mf)$weight)
  
  # Node parameters:
  # Estimate average degrees
  deg_mm <- mean(degree(net_mm))
  deg_ff <- mean(degree(net_ff))
  deg_mf <- mean(degree(net_mf)[V(net_mf)$sex == "M"])
  deg_fm <- mean(degree(net_mf)[V(net_mf)$sex == "F"])
  
  # Estimate strengths
  str_mm <- mean(strength(net_mm))
  str_ff <- mean(strength(net_ff))
  str_mf <- mean(strength(net_mf)[V(net_mf)$sex == "M"])
  str_fm <- mean(strength(net_mf)[V(net_mf)$sex == "F"])
  
  # Return a data frame
  df <- data.frame(n_edges_mm, n_edges_ff, n_edges_mf,
                   w_edges_mm, w_edges_ff, w_edges_mf,
                   deg_mm, deg_ff, deg_mf, deg_fm,
                   str_mm, str_ff, str_mf, str_fm)
  return(df)
  
}


# Observed parameters
obs_param <- getParam(net)

# Null model
ran_param <- lapply(null_net_list, getParam)
ran_param <- data.frame(rbindlist(ran_param))
ran_param


# Estimate p-values
p_val <- numeric()
for (i in 1:ncol(obs_param)) {
  if (obs_param[1, i] > mean(ran_param[, i])) {
    p <- sum(ran_param[, i] > obs_param[1, i]) / nrow(ran_param)
  } else {
    p <- sum(ran_param[, i] < obs_param[1, i]) / nrow(ran_param)
  }
  p_val <- c(p_val, p)
}
names(p_val) <- colnames(obs_param)

p_val[p_val < 0.05]
p_val[p_val > 0.05]


# Plot distributions

# Network level parameters
par(mfrow = c(2, 3))
v <- "n_edges_mm"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "n_edges_ff"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "n_edges_mf"
hist(ran_param[, v], col = "purple", border = "transparent",
     main = "M <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "w_edges_mm"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "w_edges_ff"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "w_edges_mf"
hist(ran_param[, v], col = "purple", border = "transparent",
     main = "M <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))


# Node level parameters
par(mfrow = c(2, 4))

v <- "deg_mm"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "deg_ff"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "deg_mf"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "deg_fm"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))



v <- "str_mm"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "str_ff"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "str_mf"
hist(ran_param[, v], col = "dodgerblue", border = "transparent",
     main = "M <-> F", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

v <- "str_fm"
hist(ran_param[, v], col = "coral", border = "transparent",
     main = "F <-> M", xlim = range(ran_param[, v], obs_param[, v]))
abline(v = obs_param[, v], col = ifelse(p_val[v] < 0.05, "red", "darkblue"))

