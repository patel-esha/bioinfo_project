#Hannah Luft Assignment 4

#file paths
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
diff_expr_rslts_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")

#libraries
#install.packages("mlr3")
#install.packages("e1071")
#install.packages("mlr3learners")
#install.packages("caTools") 
#install.packages("class") 
library(caTools) 
library(class) 
library(e1071)
library(magrittr)
library(dplyr)
library(mlr3)
library(mlr3learners)
#library(mlr)
#install.packages("mlr")

#load in the expression data and metadata
metadata <- readr::read_tsv(metadata_file)
# Read in the gene expression table
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Symbol")
#get rid of ham/tsp group to make it a two group problem
culledMeta <- metadata[!(metadata$refinebio_disease=="ham/tsp"),]
discardColumns <- metadata[(metadata$refinebio_disease=="ham/tsp"),]
discardColumns = as.vector(discardColumns$refinebio_accession_code)
metadata <- culledMeta
length(discardColumns)
#Preserve only columns in expression_df that match one of the accession ids
culled_expression_df = expression_df[,!(names(expression_df) %in% discardColumns)]
#check samples match (got rid of ham/tsp people)- 62 metadata cols, 63 expression (1 col for gene names)
all.equal(colnames(culled_expression_df), metadata$refinebio_accession_code)

# log scale
#epsilon <- 1e-6
#culled_expression_df <- log2(culled_expression_df + epsilon)
# order genes by variance
culled_expression_df$variance <- apply(culled_expression_df, 1, var) 
exp_ordered <- culled_expression_df[order(culled_expression_df$variance, decreasing = TRUE), ]
expressions <- select(exp_ordered, -variance)
#get top 5000
exp_top5000 <- expressions[1:5000, ]

#Supervised Analysis: 
#flip so rows are samples
flipped <- data.frame(t(exp_top5000[]))
#add disease as a column and turn them into numbers for the knn algo
flipped$disease <- metadata$refinebio_disease
flipped$disease[flipped$disease=="hc"] <- 1
flipped$disease[flipped$disease=="ms"] <- 2

#KNN

#scale
data_1 <- flipped
data_1[,1:5000] <- scale(data_1[,1:5000])
set.seed(123)

size <- floor(0.4 * nrow(data_1))
train_ind <- sample(seq_len(nrow(data_1)), size=size)
train_labels <- data_1[train_ind, 5001]
test_labels <- data_1[-train_ind, 5001]
data_train <- data_1[train_ind, 1:5000]
data_test <- data_1[-train_ind, 1:5000]
print(round(sqrt(nrow(data_train))))

predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
print(predictions)
cm <- table(test_labels, predictions)
print(cm)

#repeat with cluster numbers- grab the cluster assignments from consensusClusterPlus (saves vars in env)
combined_metadata_knn <- metadata
combined_metadata_knn$cluster_2k <- clusterAssignmentsk2g1980[[1]]
combined_metadata_knn$cluster_3k <-clusterAssignmentsk3g1980[[1]]
flipped_cluster2k <- data.frame(t(exp_top5000[]))
#add cluster as a column
flipped_cluster2k$cluster2k <- combined_metadata_knn$cluster_2k
#split data
flipped_cluster2k[,1:5000] <- scale(flipped_cluster2k[,1:5000])
size <- floor(0.4 * nrow(flipped_cluster2k))
train_ind <- sample(seq_len(nrow(flipped_cluster2k)), size=size)
train_labels <- flipped_cluster2k[train_ind, 5001]
test_labels <- flipped_cluster2k[-train_ind, 5001]
data_train <- flipped_cluster2k[train_ind, 1:5000]
data_test <- flipped_cluster2k[-train_ind, 1:5000]
print(round(sqrt(nrow(data_train))))
#do the predictions
predictions_2k <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
print(predictions_2k)
cm_2k <- table(test_labels, predictions_2k)
print(cm_2k)

#repeat with k=3
flipped_cluster3k <- data.frame(t(exp_top5000[]))
#add cluster as a column
flipped_cluster3k$cluster3k <- combined_metadata_knn$cluster_3k
#split data
flipped_cluster3k[,1:5000] <- scale(flipped_cluster3k[,1:5000])
size <- floor(0.4 * nrow(flipped_cluster3k))
train_ind <- sample(seq_len(nrow(flipped_cluster3k)), size=size)
train_labels <- flipped_cluster3k[train_ind, 5001]
test_labels <- flipped_cluster3k[-train_ind, 5001]
data_train <- flipped_cluster3k[train_ind, 1:5000]
data_test <- flipped_cluster3k[-train_ind, 1:5000]
print(round(sqrt(nrow(data_train))))
#do the predictions
predictions_3k <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
print(predictions_3k)
cm_3k <- table(test_labels, predictions_3k)
print(cm_3k)

#calculate sample area under ROC curve (AUC) across the predictive models 
#   created by each student on your team. To do this, generate a matrix of 
#   samples by models, where each cell is what class label the given model 
#   assigned to the given sample. For each sample,
#     a. How many models predict each class label, for that sample?
#     b. How many models predict the same cluster, for that sample?
#     c. Does the stability of the cluster and class label prediction 
#         correlate (calculate a statistical test to evaluate this. Don’t 
#         forget multiple test correction!)

knn_calculator_disease <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$disease <- metadata$refinebio_disease
  flipped$disease[flipped$disease=="hc"] <- 1
  flipped$disease[flipped$disease=="ms"] <- 2
  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  set.seed(123)
  #separate data
  size <- floor(0.4 * nrow(flipped))
  train_ind <- sample(seq_len(nrow(flipped)), size=size)
  train_labels <- flipped[train_ind, ngenes+1]
  test_labels <- flipped[-train_ind, ngenes+1]
  data_train <- flipped[train_ind, 1:ngenes]
  data_test <- flipped[-train_ind, 1:ngenes]
  #should be odd
  print(round(sqrt(nrow(data_train))))
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  #print(predictions)
  cm <- table(test_labels, predictions)
  print(cm)
  return(predictions)
}

knn_calculator_2k <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$cluster2k <- combined_metadata_knn$cluster_2k

  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  set.seed(123)
  #separate data
  size <- floor(0.4 * nrow(flipped))
  train_ind <- sample(seq_len(nrow(flipped)), size=size)
  train_labels <- flipped[train_ind, ngenes+1]
  test_labels <- flipped[-train_ind, ngenes+1]
  data_train <- flipped[train_ind, 1:ngenes]
  data_test <- flipped[-train_ind, 1:ngenes]
  #should be odd
  print(round(sqrt(nrow(data_train))))
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  #print(predictions)
  cm <- table(test_labels, predictions)
  print(cm)
  return(predictions)
}

knn_calculator_3k <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$cluster3k <- combined_metadata_knn$cluster_3k
  
  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  set.seed(123)
  #separate data
  size <- floor(0.4 * nrow(flipped))
  train_ind <- sample(seq_len(nrow(flipped)), size=size)
  train_labels <- flipped[train_ind, ngenes+1]
  test_labels <- flipped[-train_ind, ngenes+1]
  data_train <- flipped[train_ind, 1:ngenes]
  data_test <- flipped[-train_ind, 1:ngenes]
  #should be odd
  print(round(sqrt(nrow(data_train))))
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  #print(predictions)
  cm <- table(test_labels, predictions)
  print(cm)
  return(predictions)
}

#repeat with dif num of genes
exp_top10 <- expressions[1:10, ]
t10_disease_knn <- knn_calculator_disease(exp_top10,10)
#print(t10_disease_knn)
t10_2k_knn <- knn_calculator_2k(exp_top10, 10)
t10_3k_knn <- knn_calculator_3k(exp_top10, 10)

exp_top100 <- expressions[1:100, ]
t100_disease_knn <- knn_calculator_disease(exp_top100, 100)
#print(t100_disease_knn)
t100_2k_knn <- knn_calculator_2k(exp_top100, 100)
t100_3k_knn <- knn_calculator_3k(exp_top100, 100)

exp_top1000 <- expressions[1:1000, ]
t1000_disease_knn <- knn_calculator_disease(exp_top1000, 1000)
#print(t1000_disease_knn)
t1000_2k_knn <- knn_calculator_2k(exp_top1000, 1000)
t1000_3k_knn <- knn_calculator_3k(exp_top1000, 1000)

#AUC for each


#heatmaps and dendrograms


