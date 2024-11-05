# Define the file paths 
getwd()
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
metadata <- readr::read_tsv(metadata_file)
60


# install.packages("matrixStats")
# install.packages("factoextra")
# install.packages("ggalluvial")
# install.packages("pheatmap")
# install.packages('e1071')
# install.packages("pROC")
library(magrittr)
library(matrixStats)
library(cluster)
library("factoextra")
library(dplyr)
library(ggalluvial)
library(e1071)
library(pROC)



# Read in the gene expression table
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Symbol")
  
#get rid of ham/tsp group to make it a two group problem
culledMeta <- metadata[!(metadata$refinebio_disease=="ham/tsp"),]
discardColumns <- metadata[(metadata$refinebio_disease=="ham/tsp"),]
discardColumns <- as.vector(discardColumns$refinebio_accession_code)
metadata <- culledMeta
length(discardColumns)
#Preserve only columns in expression_df that match one of the accession ids
culled_expression_df <- expression_df[,!(names(expression_df) %in% discardColumns)]
#check samples match (got rid of ham/tsp people)- 62 metadata cols, 63 expression (1 col for gene names)
all.equal(colnames(culled_expression_df), metadata$refinebio_accession_code)

# log scale
epsilon <- 1e-6
culled_expression_df <- log2(culled_expression_df + epsilon)
# order genes by variance
culled_expression_df$variance <- apply(culled_expression_df, 1, var) 
exp_ordered <- culled_expression_df[order(culled_expression_df$variance, decreasing = TRUE), ]
expressions <- select(exp_ordered, -variance)

# Select the top 5000 rows with the highest variance
# we have to transpose since the pam method wants each row to be an observation.

exp_top10 <- expressions[1:10, ]
exp_top100 <- expressions[1:100, ]
exp_top1000 <- expressions[1:1000, ]
exp_top5000 <- expressions[1:5000, ]
exp_top10000 <- expressions[1:10000, ]
################################################ SVM for 5000 Genes
# transposes the data
transposed_df_5000 <- t(exp_top5000)
transposed_df_5000 <- as.data.frame(transposed_df_5000)

# adds a new column containing the labels (ms, hc).
labels <- metadata$refinebio_disease[match(rownames(transposed_df_5000), metadata$refinebio_accession_code)]
transposed_df_5000$Labels <- labels

# turns Labels column into categorical data
transposed_df_5000$Labels = as.factor(transposed_df_5000$Labels) 
set.seed(123) # Set some seed for reproducibility
# generating a random sample of rows for creating train dataset. Set 80% of n rows.
sample_size <- floor(0.5 * nrow(transposed_df_5000))
# creating indices for training data
train_ind <- sample(seq_len(nrow(transposed_df_5000)), size = sample_size)
# Creating the training data
train_df <- transposed_df_5000[train_ind, ]
# Creating the testing data
test_df <- transposed_df_5000[-train_ind, ]

svm_model <- svm(Labels ~ ., data = train_df, kernel = "polynomial",  cost = 1, scale = TRUE)
svm_predict_5000 <- predict(svm_model, test_df, type="response", decision.values = TRUE)

svm_predict_vs_true_5000 <- data.frame(test_df$Labels, svm_predict_5000)

# Convert to a data frame suitable for ggplot2
names(svm_predict_vs_true_5000) <- c("True", "Predicted")

# Plot using ggplot2
png("plots/assignment4/svm_5000_ggplot2.png")
ggplot(data = svm_predict_vs_true_5000, aes(x = True, y = Predicted, color = True)) +
    geom_jitter(width = 0.2, height = 0.2, shape = 16, alpha = 0.6) +
    labs(title = "SVM Prediction with 5000 Genes", x = "True Categories", y = "Predicted Categories") +
    theme_minimal() +
    scale_color_manual(values = c("ms" = "blue", "hc" = "red"))
dev.off()

print(mean(svm_predict_5000 == test_df$Labels))

label_mapping <- c("ms"="1", "hc"="2")
svm_predict_vs_true_5000$True <- as.numeric(label_mapping[as.character(svm_predict_vs_true_5000$True)])
svm_predict_vs_true_5000$Predicted <- as.numeric(label_mapping[as.character(svm_predict_vs_true_5000$Predicted)])
roc_curve <- roc(svm_predict_vs_true_5000$True, svm_predict_vs_true_5000$Predicted)
auc_value <- auc(roc_curve)
auc_value

########################################## SVM for clusters from Assignment 3


# transposes the data
transposed_df_5000 <- t(exp_top5000)
transposed_df_5000 <- as.data.frame(transposed_df_5000)
label_mapping <- c("1" = "ms", "2" = "hc")
labels_df <- data.frame(
  refinebio_accession_code = names(assn3ClusterResult),
  Labels = label_mapping[as.character(assn3ClusterResult)],
  stringsAsFactors = FALSE
)


# adds a new column containing the labels (ms, hc).
labels <- metadata$refinebio_disease[match(rownames(transposed_df_5000), labels_df$refinebio_accession_code)]
transposed_df_5000$Labels <- labels

# turns Labels column into categorical data
transposed_df_5000$Labels = as.factor(transposed_df_5000$Labels) 
set.seed(123) # Set some seed for reproducibility
# generating a random sample of rows for creating train dataset. Set 80% of n rows.
sample_size <- floor(0.8 * nrow(transposed_df_5000))
# creating indices for training data
train_ind <- sample(seq_len(nrow(transposed_df_5000)), size = sample_size)
# Creating the training data
train_df <- transposed_df_5000[train_ind, ]
# Creating the testing data
test_df <- transposed_df_5000[-train_ind, ]

svm_model_assn3 <- svm(Labels ~ ., data = train_df, kernel = "linear",  cost = 1, scale = TRUE)
svm_predict_5000_assn3 <- predict(svm_model_assn3, test_df, type="response", decision.values = TRUE)

# Plot using ggplot2
png("plots/assignment4/svm_5000_assn3_clust_ggplot2.png")
ggplot(data = svm_predict_vs_true_5000, aes(x = True, y = Predicted, color = True)) +
  geom_jitter(width = 0.2, height = 0.2, shape = 16, alpha = 0.6) +
  labs(title = "SVM Prediction for clusters from Assignment 3 with 5000 Genes", x = "True Categories", y = "Predicted Categories") +
  theme_minimal() +
  scale_color_manual(values = c("ms" = "blue", "hc" = "red"))
dev.off()

print(mean(svm_predict_5000_assn3 == test_df$Labels))


################################################ SVM for 10 Genes
# transposes the data
transposed_df_10 <- t(exp_top10)
transposed_df_10 <- as.data.frame(transposed_df_10)

# adds a new column containing the labels (ms, hc).
labels <- metadata$refinebio_disease[match(rownames(transposed_df_10), metadata$refinebio_accession_code)]
transposed_df_10$Labels <- labels

# turns Labels column into categorical data
transposed_df_10$Labels = as.factor(transposed_df_10$Labels) 
set.seed(123) # Set some seed for reproducibility
# generating a random sample of rows for creating train dataset. Set 80% of n rows.
sample_size <- floor(0.6 * nrow(transposed_df_10))
# creating indices for training data
train_ind <- sample(seq_len(nrow(transposed_df_10)), size = sample_size)
# Creating the training data
train_df <- transposed_df_10[train_ind, ]
# Creating the testing data
test_df <- transposed_df_10[-train_ind, ]

svm_model_10 <- svm(Labels ~ ., data = train_df, kernel = "polynomial",  cost = 1, scale = TRUE)
svm_predict_10 <- predict(svm_model_10, test_df, type="response", decision.values = TRUE)

svm_predict_vs_true_10 <- data.frame(test_df$Labels, svm_predict_10)

# Convert to a data frame suitable for ggplot2
names(svm_predict_vs_true_10) <- c("True", "Predicted")

# Plot using ggplot2
png("plots/assignment4/svm_10_ggplot2.png")
ggplot(data = svm_predict_vs_true_10, aes(x = True, y = Predicted, color = True)) +
  geom_jitter(width = 0.2, height = 0.2, shape = 16, alpha = 0.6) +
  labs(title = "SVM Prediction with 10 Genes", x = "True Categories", y = "Predicted Categories") +
  theme_minimal() +
  scale_color_manual(values = c("ms" = "blue", "hc" = "red"))
dev.off()

print(mean(svm_predict_10 == test_df$Labels))

label_mapping <- c("ms"="1", "hc"="2")
svm_predict_vs_true_10$True <- as.numeric(label_mapping[as.character(svm_predict_vs_true_10$True)])
svm_predict_vs_true_10$Predicted <- as.numeric(label_mapping[as.character(svm_predict_vs_true_10$Predicted)])
roc_curve <- roc(svm_predict_vs_true_10$True, svm_predict_vs_true_10$Predicted)
auc_value <- auc(roc_curve)
auc_value
svm_model_10

################################################ SVM for 100 Genes
# transposes the data
transposed_df_100 <- t(exp_top100)
transposed_df_100 <- as.data.frame(transposed_df_100)

# adds a new column containing the labels (ms, hc).
labels <- metadata$refinebio_disease[match(rownames(transposed_df_100), metadata$refinebio_accession_code)]
transposed_df_100$Labels <- labels

# turns Labels column into categorical data
transposed_df_100$Labels = as.factor(transposed_df_100$Labels) 
set.seed(123) # Set some seed for reproducibility
# generating a random sample of rows for creating train dataset. Set 80% of n rows.
sample_size <- floor(0.5 * nrow(transposed_df_100))
# creating indices for training data
train_ind <- sample(seq_len(nrow(transposed_df_100)), size = sample_size)
# Creating the training data
train_df <- transposed_df_100[train_ind, ]
# Creating the testing data
test_df <- transposed_df_100[-train_ind, ]

svm_model_100 <- svm(Labels ~ ., data = train_df, kernel = "polynomial",  cost = 1, scale = TRUE)
svm_predict_100 <- predict(svm_model, test_df, type="response", decision.values = TRUE)

svm_predict_vs_true_100 <- data.frame(test_df$Labels, svm_predict_100)

# Convert to a data frame suitable for ggplot2
names(svm_predict_vs_true_100) <- c("True", "Predicted")

# Plot using ggplot2
png("plots/assignment4/svm_100_ggplot2.png")
ggplot(data = svm_predict_vs_true_100, aes(x = True, y = Predicted, color = True)) +
  geom_jitter(width = 0.2, height = 0.2, shape = 16, alpha = 0.6) +
  labs(title = "SVM Prediction with 100 Genes", x = "True Categories", y = "Predicted Categories") +
  theme_minimal() +
  scale_color_manual(values = c("ms" = "blue", "hc" = "red"))
dev.off()

print(mean(svm_predict_100 == test_df$Labels))

label_mapping <- c("ms"="1", "hc"="2")
svm_predict_vs_true_100$True <- as.numeric(label_mapping[as.character(svm_predict_vs_true_100$True)])
svm_predict_vs_true_100$Predicted <- as.numeric(label_mapping[as.character(svm_predict_vs_true_100$Predicted)])
roc_curve <- roc(svm_predict_vs_true_100$True, svm_predict_vs_true_100$Predicted)
auc_value <- auc(roc_curve)
auc_value
svm_model_100
################################################ SVM for 10000 Genes
# transposes the data
transposed_df_10000 <- t(exp_top10000)
transposed_df_10000 <- as.data.frame(transposed_df_10000)

# adds a new column containing the labels (ms, hc).
labels <- metadata$refinebio_disease[match(rownames(transposed_df_10000), metadata$refinebio_accession_code)]
transposed_df_10000$Labels <- labels

# turns Labels column into categorical data
transposed_df_10000$Labels = as.factor(transposed_df_10000$Labels) 
set.seed(123) # Set some seed for reproducibility
# generating a random sample of rows for creating train dataset. Set 80% of n rows.
sample_size <- floor(0.6 * nrow(transposed_df_10000))
# creating indices for training data
train_ind <- sample(seq_len(nrow(transposed_df_10000)), size = sample_size)
# Creating the training data
train_df <- transposed_df_10000[train_ind, ]
# Creating the testing data
test_df <- transposed_df_10000[-train_ind, ]

rownames(test_df)

svm_model_10000 <- svm(Labels ~ ., data = train_df, kernel = "polynomial",  cost = 1, scale = TRUE)
svm_predict_10000 <- predict(svm_model_10000, test_df, type="response", decision.values = TRUE)

svm_predict_vs_true_10000 <- data.frame(test_df$Labels, svm_predict_10000)

# Convert to a data frame suitable for ggplot2
names(svm_predict_vs_true_10000) <- c("True", "Predicted")

# Plot using ggplot2
png("plots/assignment4/svm_10000_ggplot2.png")
ggplot(data = svm_predict_vs_true_10000, aes(x = True, y = Predicted, color = True)) +
  geom_jitter(width = 0.2, height = 0.2, shape = 16, alpha = 0.6) +
  labs(title = "SVM Prediction with 10000 Genes", x = "True Categories", y = "Predicted Categories") +
  theme_minimal() +
  scale_color_manual(values = c("ms" = "blue", "hc" = "red"))
dev.off()

print(mean(svm_predict_10000 == test_df$Labels))

label_mapping <- c("ms"="1", "hc"="2")
svm_predict_vs_true_10000$True <- as.numeric(label_mapping[as.character(svm_predict_vs_true_10000$True)])
svm_predict_vs_true_10000$Predicted <- as.numeric(label_mapping[as.character(svm_predict_vs_true_10000$Predicted)])
roc_curve <- roc(svm_predict_vs_true_10000$True, svm_predict_vs_true_10000$Predicted)
auc_value <- auc(roc_curve)
auc_value
svm_model_10000

################################ Heatmap with dendogram
culledMeta <- metadata[!(metadata$refinebio_disease=="ham/tsp"),]
discardColumns <- metadata[(metadata$refinebio_disease=="ham/tsp"),]
discardColumns = as.vector(discardColumns$refinebio_accession_code)
metadata <- culledMeta
length(discardColumns)
coefs <- t(svm_model_10$coefs) %*% svm_model_10$SV
feature_importance <- data.frame(t(coefs))
feature_importance <- feature_importance %>%
  filter(abs(t.coefs.) >= 10)
# filters the impactful genes from the data set
common_row_names_svm <- intersect(rownames(feature_importance), rownames(exp_top10))
filtered_expression_df_svm <- expression_df[common_row_names_svm, ]

# reads in meta data with rowname already set
metadata_withRowNames <- metadata %>%
  tibble::column_to_rownames("refinebio_accession_code")



# creates the annotation data frame
annotation <- metadata_withRowNames[c("refinebio_disease")]
annotation$refinebio_disease <- ifelse(grepl("^hc", annotation$refinebio_disease), "hc",
                                            ifelse(grepl("^ms ", annotation$refinebio_disease), "ms", annotation$refinebio_disease))
names(annotation)[names(annotation) == "refinebio_diesease"] <- "SampleGroup"

# Heatmap of svm with 100 genes

png("plots/assignment4/heatmap_svm_100genes.png")

# Create the heatmap
heatmap <- pheatmap(
  filtered_expression_df_svm,
  annotation_col = annotation,
  show_rownames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "none",
  main = "Impactfull Genes",
)

print(heatmap)
dev.off()

