# Assignment 4- Peter Walsh
library(magrittr)
library(dplyr)
library(tidymodels)
library(readr)



data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
diff_expr_rslts_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")



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




culled_expression_df$variance <- apply(culled_expression_df, 1, var) 
exp_ordered <- culled_expression_df[order(culled_expression_df$variance, decreasing = TRUE), ]
expressions <- select(exp_ordered, -variance)
exp_top5000 <- expressions[1:5000, ]


# Prepare Data
# Define the target variable
metadata$group <- ifelse(metadata$refinebio_disease == "hc", 0, 1)
# Make sure both datasets align in order
metadata <- metadata[match(colnames(exp_top5000), metadata$refinebio_accession_code), ]
all.equal(colnames(exp_top5000), metadata$refinebio_accession_code)



# Split into training and testing sets
# Combine expression data with the group label
exp_data <- t(exp_top5000)
exp_data <- as.data.frame(exp_data)
exp_data$group <- metadata$group

# Create train-test split
set.seed(123)  # For reproducibility
data_split <- initial_split(exp_data, prop = 0.8, strata = group)
train_data <- training(data_split)
test_data <- testing(data_split)
# Convert 'group' to a factor in both train and test sets
train_data$group <- as.factor(train_data$group)
test_data$group <- as.factor(test_data$group)


# Define and Train Logistic Regression Model
# Define the logistic regression model
log_reg_model <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")
# Define a workflow
log_reg_workflow <- workflow() %>%
  add_model(log_reg_model) %>%
  add_formula(group ~ .)
# Fit the model
log_reg_fit <- log_reg_workflow %>%
  fit(data = train_data)




# Evaluate Model Performance on test set
# Predict on the test data
test_predictions <- predict(log_reg_fit, test_data) %>%
  bind_cols(test_data) %>%
  metrics(truth = group, estimate = .pred_class)
# View test set performance
print(test_predictions)



# Retrain Model for Multi-Class Cluster Prediction
metadata$cluster <- clusterAssignmentsk2g1980[[1]]

# Update training data to use clusters as the target
exp_data$cluster <- metadata$cluster

# Split data for multi-class prediction
data_split <- initial_split(exp_data, prop = 0.8, strata = cluster)
train_data <- training(data_split)
test_data <- testing(data_split)
# Convert 'cluster' to a factor in both train and test sets
train_data$cluster <- as.factor(train_data$cluster)
exp_data$cluster <- as.factor(exp_data$cluster)


# Update workflow for multi-class prediction
log_reg_workflow_multi <- workflow() %>%
  add_model(log_reg_model) %>%
  add_formula(cluster ~ .)

# Fit model to multi-class data
log_reg_fit_multi <- log_reg_workflow_multi %>%
  fit(data = train_data)

# Predict and evaluate on the test set
test_predictions_multi <- predict(log_reg_fit_multi, test_data) %>%
  bind_cols(test_data) %>%
  metrics(truth = cluster, estimate = .pred_class)

print(test_predictions_multi)

predictions <- predict(log_reg_fit, exp_data) %>%
  bind_cols(exp_data)




