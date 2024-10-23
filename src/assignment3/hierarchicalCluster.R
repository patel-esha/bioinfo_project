# Hierarchical Clustering- Peter Walsh
library(cluster)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)

data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
diff_expr_rslts_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")

# Load in the expression data and metadata
expressions <- readr::read_tsv(data_file)
metadata <- readr::read_tsv(metadata_file)
diff_expr <- readr::read_tsv(diff_expr_rslts_file)

# Extend expressions with their variance
expressions$variance <- apply(expressions[,-1], 1, var) 
exp_ordered <- expressions[order(-expressions$variance), ]

# Function to perform hierarchical clustering and silhouette analysis
optimal_clusters_analysis <- function(data, max_clusters) {
  data_label <- data$Symbol
  data$Symbol <- NULL
  data_sc <- as.data.frame(scale(data))
  dist_mat <- dist(data_sc, method='euclidean')
  hclust_res <- hclust(dist_mat, method='average')
  
  silhouette_scores <- numeric(max_clusters)
  cluster_assignments <- list()  # Store cluster assignments for chi-squared tests
  
  for (k in 2:max_clusters) {
    clusters <- cutree(hclust_res, k = k)
    
    # Store cluster assignments
    cluster_assignments[[k]] <- clusters
    
    # Ensure there are at least two clusters
    if (length(unique(clusters)) > 1) {
      sil <- silhouette(clusters, dist_mat)
      if (is.matrix(sil)) {
        silhouette_scores[k] <- mean(sil[, 3], na.rm = TRUE)  # Average silhouette width
      } else {
        silhouette_scores[k] <- NA  # If sil is not a matrix, assign NA
      }
    } else {
      silhouette_scores[k] <- NA  # Not enough clusters to compute silhouette
    }
  }
  
  # Handle NA values before finding the optimal number of clusters
  valid_scores <- silhouette_scores[!is.na(silhouette_scores)]
  optimal_k <- if (length(valid_scores) > 0) {
    which.max(silhouette_scores)  # Find the index of the max value
  } else {
    NA  # If no valid scores, return NA
  }
  
  return(list(optimal_k = optimal_k, silhouette_scores = silhouette_scores, cluster_assignments = cluster_assignments))
}

# List of top gene counts to analyze
top_gene_counts <- c(10, 100, 1000, 5000, 10000)

# Store optimal clusters, silhouette scores, and assignments for each top gene count
results <- lapply(top_gene_counts, function(n) {
  optimal_clusters_analysis(exp_ordered[1:n, ], max_clusters = 10)  # Adjust max_clusters as needed
})

# Combine results into a data frame for easier visualization
results_df <- data.frame(
  Top_Genes = top_gene_counts,
  Optimal_Clusters = sapply(results, function(x) x$optimal_k),
  Silhouette_Scores = sapply(results, function(x) max(x$silhouette_scores, na.rm = TRUE))
)

print(results_df)


