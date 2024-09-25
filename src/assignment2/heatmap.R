#heatmap

# Define the file paths 
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
top_50_file <- file.path(results_dir, "SRP164913_diff_expr_top_50_results.tsv")


#install once
#install_github("jokergoo/ComplexHeatmap")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#libs
library(devtools)
library(ComplexHeatmap)
library(circlize)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(dplyr)


#read in the top 50
top_50 <- readr::read_tsv(top_50_file)
#read in the hugo data and grab data from the top 50 genes? 
expression_df <- readr::read_tsv(data_file)

#read in metadata
metadata <- readr::read_tsv(metadata_file)

#get rid of ham/tsp
culledMeta <- metadata[!(metadata$refinebio_disease=="ham/tsp"),]
discardColumns <- metadata[(metadata$refinebio_disease=="ham/tsp"),]
discardColumns = as.vector(discardColumns$refinebio_accession_code)
metadata <- culledMeta
length(discardColumns)
#Preserve only columns in expression_df that match one of the accession ids
culled_expression_df = expression_df[,!(names(expression_df) %in% discardColumns)]
#check samples match (got rid of ham/tsp people)
all.equal(colnames(culled_expression_df), metadata$refinebio_accession_code)

#grab only top 50 genes
gene_list <- as.list(top_50[1])
#print(gene_list)
#print(expression_df)
list_i <- c()
for(i in 1:nrow(expression_df)){
  if(expression_df[i,1] %in% top_50$Gene){
    list_i <- c(list_i, i)
    #print(TRUE)
  }
}
top_50_expression_df <- culled_expression_df[list_i,]

metadata <- metadata %>%
  dplyr::mutate(
    refinebio_disease = factor(
      refinebio_disease,
      # specify the possible levels in the order we want them to appear
      levels = c("hc", "ms")
    )
  )

#create heatmap
#remove symbol column
top_50_expression_df$Symbol <- NULL
expression_matrix <- data.matrix(top_50_expression_df)

#"add a sidebar colored by sample groupings"
column_ha = HeatmapAnnotation(labels = as.factor(metadata$refinebio_disease))

png("SRP164913_heatmap_plot.png")
ht <- Heatmap(
  expression_matrix, 
  name="Top 50", 
  top_annotation = column_ha,
  column_title = "Samples", column_title_side = "bottom", 
  row_title = "Genes", row_title_side = "right" ,
  show_column_names = FALSE
)
draw(ht)
dev.off()
