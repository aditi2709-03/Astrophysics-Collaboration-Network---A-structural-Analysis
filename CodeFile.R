library(igraph)
library(tidyverse)

astro_ph <- read_graph("astro-ph.gml", format = "gml")

plot(astro_ph, vertex.size = 5, vertex.label = NA, edge.arrow.size = 0.5, 
     main = "Astrophysics Collaboration Network")

# Create a smaller subgraph for visualization
set.seed(59)
sample_nodes <- sample(V(astro_ph), size = min(500, vcount(astro_ph)))  # Select 500 nodes or fewer
subgraph <- induced_subgraph(astro_ph, vids = sample_nodes) 
plot(subgraph, vertex.size = 2, vertex.label = NA, edge.width = 0.5,edge.color = "black")

vcount(astro_ph) #nodes
ecount(astro_ph) #edges
mean(degree(astro_ph)) #avg degree 

# Compute degree of each node -- Degree Distribution 
node_degrees <- degree(astro_ph)
# Create degree distribution
degree_dist <- table(node_degrees) / vcount(astro_ph)
degree_df <- data.frame(k = as.numeric(names(degree_dist)), Pk = as.numeric(degree_dist))

# Plot degree distribution (log-log scale)
log_log <- ggplot(degree_df, aes(x = k, y = Pk)) +
  geom_point(color = "blue") +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Degree Distribution (Log-Log Plot)", x = "Degree (k)", y = "P(k)") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("log-log.png", plot = log_log, width = 9, height = 4, dpi = 300)

# Shortest path and mean distance
shortest_paths <- distances(astro_ph)
mean_distance <- mean(shortest_paths[shortest_paths != Inf])
cat("Mean shortest path length:", mean_distance, "\n")

# Assuming your network is stored in a variable called 'g'
components_info <- components(astro_ph)

# Number of connected components
num_components <- components_info$no
cat("Number of connected components:", num_components, "\n")

# Sizes of each component
component_sizes <- components_info$csize
print(component_sizes)

# Find the largest component (the one with the maximum size)
largest_comp_id <- which.max(components_info$csize)

# Extract nodes belonging to the largest component
largest_comp_nodes <- which(components_info$membership == largest_comp_id)

# Induce subgraph with only these nodes
largest_subgraph <- induced_subgraph(astro_ph, vids = largest_comp_nodes)

# Plot the largest component
plot(largest_subgraph, 
     vertex.size = 5, 
     vertex.label = NA,  # Remove labels for clarity
     edge.arrow.size = 0.5, 
     main = "Largest Connected Component")


# --- CLUSTERING COEFFICIENT ANALYSIS ---

# Compute global clustering coefficient
global_clustering <- transitivity(astro_ph, type = "global")

# Compute local clustering coefficient for each node
local_clustering <- transitivity(astro_ph, type = "local", isolates = "zero")

# Summary statistics for clustering coefficients
clustering_summary <- tibble(node = names(V(astro_ph)), clustering = local_clustering) %>%
  summarise(mean_clustering = mean(clustering, na.rm = TRUE),
            median_clustering = median(clustering, na.rm = TRUE))

print(paste("Global Clustering Coefficient:", round(global_clustering, 4)))
print(clustering_summary)

# Plot Distribution of Local Clustering Coefficients
CC <- ggplot(tibble(clustering = local_clustering), aes(x = clustering)) +
  geom_histogram(bins = 10, fill = "grey", color = "black") +
  xlab("Local Clustering Coefficient") +
  ylab("Frequency") +
  ggtitle("Distribution of Local Clustering Coefficients") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("CC.png", plot = CC, width = 8, height = 4, dpi = 300)



# rp ----------------------------------------------------------------------

# --- CENTRALITY ANALYSIS ---

# Degree centrality (already implicit, but make it explicit)
degree_centrality <- degree(astro_ph)

# Betweenness centrality – identifies boundary spanners
betweenness_centrality <- betweenness(astro_ph, directed = FALSE)

# Closeness centrality – nodes that can quickly interact
closeness_centrality <- closeness(astro_ph)

# Eigenvector centrality – influence based on influential neighbors
eigen_centrality <- eigen_centrality(astro_ph)$vector

# Combine into a centrality dataframe
centrality_df <- tibble(
  node = names(V(astro_ph)),
  degree = degree_centrality,
  betweenness = betweenness_centrality,
  closeness = closeness_centrality,
  eigenvector = eigen_centrality
)

# Save to CSV for reference
write.csv(centrality_df, "centrality_metrics.csv", row.names = FALSE)

# --- ASSORTATIVITY ---

assortativity_val <- assortativity_degree(astro_ph, directed = FALSE)
cat("Assortativity coefficient (degree correlation):", assortativity_val, "\n")

# --- COMMUNITY DETECTION ---

louvain_clusters <- cluster_louvain(astro_ph)
cat("Number of Louvain communities:", length(louvain_clusters), "\n")

# Plot with communities
png("louvain_communities.png", width = 1200, height = 1000)
plot(louvain_clusters, astro_ph, vertex.size = 2, vertex.label = NA,
     main = "Louvain Communities in astro-ph")
dev.off()  



# --- CLUSTERING VS DEGREE ---

clust_deg_df <- tibble(
  degree = degree_centrality,
  clustering = local_clustering
)

gg_clust_deg <- ggplot(clust_deg_df, aes(x = degree, y = clustering)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  labs(title = "Clustering Coefficient vs. Degree",
       x = "Degree (log scale)",
       y = "Local Clustering Coefficient") +
  theme_minimal()
ggsave("clustering_vs_degree.png", plot = gg_clust_deg, width = 8, height = 5, dpi = 300)


# --- DEGREE DISTRIBUTION FITTING ---
install.packages("poweRlaw")
library(poweRlaw)

# Filter out 0s
degree_nonzero <- degree_centrality[degree_centrality > 0]

# Fit using poweRlaw
pl_model <- displ$new(degree_nonzero)
est <- estimate_xmin(pl_model)
pl_model$setXmin(est)

png("powerlaw_fit.png", width = 800, height = 600)
plot(pl_model)
lines(pl_model, col = "red")
dev.off()

cat("Fitted Power-law alpha:", est$pars, "\n")

