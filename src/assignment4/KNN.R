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
library(ComplexHeatmap)
#install.packages("pROC")
library(pROC)
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
set.seed(123)

#Supervised Analysis: 
#get the cluster groups from consensusClusterPlus
combined_metadata_knn <- metadata
combined_metadata_knn$cluster_2k <- clusterAssignmentsk2g1980[[1]]
combined_metadata_knn$cluster_3k <-clusterAssignmentsk3g1980[[1]]
knn_calculator_disease <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$disease <- metadata$refinebio_disease
  flipped$disease[flipped$disease=="hc"] <- 1
  flipped$disease[flipped$disease=="ms"] <- 2
  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  #separate data
  #size <- floor(0.5 * nrow(flipped))
  #train_ind <- sample(seq_len(nrow(flipped)), size=size)
  #train_labels <- flipped[train_ind, ngenes+1]
  #test_labels <- flipped[-train_ind, ngenes+1]
  #data_train <- flipped[train_ind, 1:ngenes]
  #data_test <- flipped[-train_ind, 1:ngenes]
  
  test_samples <- c("SRR7993430", "SRR7993431", "SRR7993433", "SRR7993436", "SRR7993441", "SRR7993447", "SRR7993449", "SRR7993451", "SRR7993452", "SRR7993454", "SRR7993455", "SRR7993457", "SRR7993495", "SRR7993499", "SRR7993500", "SRR7993501", "SRR7993507", "SRR7993508", "SRR7993514", "SRR7993515", "SRR7993521", "SRR7993524", "SRR7993525", "SRR7993526", "SRR7993527")
  data_test <- subset(flipped, rownames(flipped) %in% test_samples)
  data_train <- subset(flipped, !(rownames(flipped) %in% test_samples))
  train_labels <- data_train[,ngenes+1]
  test_labels <- data_test[,ngenes+1]
  data_test <- data_test[,1:ngenes]
  data_train <- data_train[,1:ngenes]
  
  
  #should be odd
  #print(round(sqrt(nrow(data_train))))
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  print(predictions)
  #predictionFull <- flipped
  #predictionFull$results <- flipped$disease
  #k=1
  #print(length(predictions))
  #for(i in nrow(predictionFull)){
  #  if(i %in% train_ind){
  #    
  #  }
  #  else{
  #    predictionFull[i, ngenes+2] <- predictions[k]
  #  }
  #  k <- k+1
  #}
  #print(predictionFull$results)
  cm <- table(test_labels, predictions)
  print(cm)
  return(list(predictions, cm, test_labels, predictions))
}

knn_calculator_2k <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$cluster2k <- combined_metadata_knn$cluster_2k
  
  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  #separate data
  #size <- floor(0.5 * nrow(flipped))
  #train_ind <- sample(seq_len(nrow(flipped)), size=size)
  #train_labels <- flipped[train_ind, ngenes+1]
  #test_labels <- flipped[-train_ind, ngenes+1]
  #data_train <- flipped[train_ind, 1:ngenes]
  #data_test <- flipped[-train_ind, 1:ngenes]
  #should be odd
  #print(round(sqrt(nrow(data_train))))
  test_samples <- c("SRR7993430", "SRR7993431", "SRR7993433", "SRR7993436", "SRR7993441", "SRR7993447", "SRR7993449", "SRR7993451", "SRR7993452", "SRR7993454", "SRR7993455", "SRR7993457", "SRR7993495", "SRR7993499", "SRR7993500", "SRR7993501", "SRR7993507", "SRR7993508", "SRR7993514", "SRR7993515", "SRR7993521", "SRR7993524", "SRR7993525", "SRR7993526", "SRR7993527")
  data_test <- subset(flipped, rownames(flipped) %in% test_samples)
  data_train <- subset(flipped, !(rownames(flipped) %in% test_samples))
  train_labels <- data_train[,ngenes+1]
  test_labels <- data_test[,ngenes+1]
  data_test <- data_test[,1:ngenes]
  data_train <- data_train[,1:ngenes]
  
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  #print(predictions)
  # predictionFull <- flipped
  # predictionFull$results <- flipped$cluster2k
  # k=1
  # print(length(predictions))
  # for(i in nrow(predictionFull)){
  #   if(i %in% train_ind){
  #     
  #   }
  #   else{
  #     predictionFull[i, ngenes+2] <- predictions[k]
  #   }
  #   k <- k+1
  # }
  #print(predictionFull$results)
  cm <- table(test_labels, predictions)
  print(cm)
  return(list(predictions, cm, test_labels, predictions))
}

knn_calculator_3k <- function(exp_df, ngenes){
  #flip the data
  flipped <- data.frame(t(exp_df[]))
  #add disease as a column and turn them into numbers for the knn algo
  flipped$cluster3k <- combined_metadata_knn$cluster_3k
  
  #scale data
  flipped[,1:ngenes] <- scale(flipped[,1:ngenes])
  #separate data
  #size <- floor(0.5 * nrow(flipped))
  #train_ind <- sample(seq_len(nrow(flipped)), size=size)
  #train_labels <- flipped[train_ind, ngenes+1]
  #test_labels <- flipped[-train_ind, ngenes+1]
  #data_train <- flipped[train_ind, 1:ngenes]
  #data_test <- flipped[-train_ind, 1:ngenes]
  
  test_samples <- c("SRR7993430", "SRR7993431", "SRR7993433", "SRR7993436", "SRR7993441", "SRR7993447", "SRR7993449", "SRR7993451", "SRR7993452", "SRR7993454", "SRR7993455", "SRR7993457", "SRR7993495", "SRR7993499", "SRR7993500", "SRR7993501", "SRR7993507", "SRR7993508", "SRR7993514", "SRR7993515", "SRR7993521", "SRR7993524", "SRR7993525", "SRR7993526", "SRR7993527")
  data_test <- subset(flipped, rownames(flipped) %in% test_samples)
  data_train <- subset(flipped, !(rownames(flipped) %in% test_samples))
  train_labels <- data_train[,ngenes+1]
  test_labels <- data_test[,ngenes+1]
  data_test <- data_test[,1:ngenes]
  data_train <- data_train[,1:ngenes]
  
  
  #should be odd
  print(round(sqrt(nrow(data_train))))
  
  #do KNN
  predictions <- knn(train = data_train, test = data_test, cl = train_labels,k = round(sqrt(nrow(data_train))))
  #print(predictions)
  # predictionFull <- flipped
  # predictionFull$results <- flipped$cluster3k
  # k=1
  # print(length(predictions))
  # for(i in nrow(predictionFull)){
  #   if(i %in% train_ind){
  #     
  #   }
  #   else{
  #     predictionFull[i, ngenes+2] <- predictions[k]
  #   }
  #   k <- k+1
  # }
  #print(predictionFull$results)
  cm <- table(test_labels, predictions)
  print(cm)
  return(list(predictions, cm, test_labels, predictions))
}


#KNN
knn_disease_og <- knn_calculator_disease(exp_top5000, 5000)

#repeat with cluster numbers
knn_2k_og <- knn_calculator_2k(exp_top5000, 5000)
knn_3k_og <- knn_calculator_3k(exp_top5000, 5000)


#calculate sample area under ROC curve (AUC) across the predictive models 
#   created by each student on your team. To do this, generate a matrix of 
#   samples by models, where each cell is what class label the given model 
#   assigned to the given sample. For each sample,
#     a. How many models predict each class label, for that sample?
#     b. How many models predict the same cluster, for that sample?
#     c. Does the stability of the cluster and class label prediction 
#         correlate (calculate a statistical test to evaluate this. Don’t 
#         forget multiple test correction!)


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

#pROC
roc_5000disease <- roc(knn_disease_og[[3]], as.numeric(knn_disease_og[[4]]))
auc1 <- auc(roc_5000disease)
auc1
roc_50002k <- roc(knn_2k_og[[3]], as.numeric(knn_2k_og[[4]]))
auc <- auc(roc_50002k)
auc
roc_50003k <- roc(knn_3k_og[[3]], as.numeric(knn_3k_og[[4]]))
auc <- auc(roc_50003k)
auc

roc_1000disease <- roc(t1000_disease_knn[[3]], as.numeric(t1000_disease_knn[[4]]))
auc <- auc(roc_1000disease)
auc
roc_10002k <- roc(t1000_2k_knn[[3]], as.numeric(t1000_2k_knn[[4]]))
auc <- auc(roc_10002k)
auc
roc_10003k <-roc(t1000_3k_knn[[3]], as.numeric(t1000_3k_knn[[4]]))
auc <- auc(roc_10003k)
auc

roc_100disease <- roc(t100_disease_knn[[3]], as.numeric(t100_disease_knn[[4]]))
auc <- auc(roc_100disease)
auc
roc_1002k <- roc(t100_2k_knn[[3]], as.numeric(t100_2k_knn[[4]]))
auc <- auc(roc_1002k)
auc
roc_1003k <- roc(t100_3k_knn[[3]], as.numeric(t100_3k_knn[[4]]))
auc <- auc(roc_1003k)
auc


roc_10disease <- roc(t10_disease_knn[[3]], as.numeric(t10_disease_knn[[4]]))
print(t10_disease_knn[[4]])
auc <- auc(roc_10disease)
auc
roc_102k <- roc(t10_2k_knn[[3]], as.numeric(t10_2k_knn[[4]]))
auc <- auc(roc_t102k)
auc
roc_103k <- roc(t10_3k_knn[[3]], as.numeric(t10_3k_knn[[4]]))
auc <- auc(roc_103k)
auc

#heatmaps and dendrograms
#put together predictions / training based on indices for the knn results 
print(knn_disease_og[3])
print(length(t(knn_disease_og[3])[1]))


column_ha = HeatmapAnnotation(groups = as.factor(metadata$refinebio_disease), predictedgroups = as.factor(t(knn_disease_og[3])[[1]]), predicted2k = as.factor(t(knn_2k_og[3])[[1]]), predicted3k = as.factor(t(knn_3k_og[3])[[1]]))

png("results/knnResults/knnHeatmap.png")

ht <- Heatmap(
  as.matrix(exp_top100), 
  name="KNN Heatmap, 100 Genes", 
  top_annotation = column_ha,
  column_title = "Samples", column_title_side = "bottom", 
  row_title = "Genes", row_title_side = "right" ,
  show_column_names = FALSE,
  show_row_names = FALSE
)

draw(ht)
dev.off()


