#gprofiler2
# Define the file paths 
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
diff_expr_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")
plots_dir <- file.path("plots")


#installs
#install.packages("gprofiler2")
#libs
library(gprofiler2)

#load in raw data
#expression_df <- readr::read_tsv(data_file)
#gene_list <- as.list(expression_df[1])
#cull data
#read in metadata
#metadata <- readr::read_tsv(metadata_file)
#get rid of ham/tsp
#culledMeta <- metadata[!(metadata$refinebio_disease=="ham/tsp"),]
#discardColumns <- metadata[(metadata$refinebio_disease=="ham/tsp"),]
#discardColumns = as.vector(discardColumns$refinebio_accession_code)
#metadata <- culledMeta
#length(discardColumns)
#Preserve only columns in expression_df that match one of the accession ids
#culled_expression_df = expression_df[,!(names(expression_df) %in% discardColumns)]
#check samples match (got rid of ham/tsp people)
#all.equal(colnames(culled_expression_df), metadata$refinebio_accession_code)


#load in differential expression data
diff_expr_df <- readr::read_tsv(diff_expr_file)

#ontology: ms vs healthy control (disease)

#enrichment analysis gene ontology
#gost_res <- gost(
#  query = gene_list[1:50],  #list of gene names (hugo)
#  organism = "hsapiens",  #human species
#  ordered_query = FALSE, #the data is not ordered
#  significant = FALSE, #want everything, not just the statistically significant
#  exclude_iea = TRUE, #true to do gene oncology
#  measure_underrepresentation = FALSE, #don't measure underrepresentation
#  evcodes = FALSE, #
#  user_threshold = 0.05, #default for the pvalue threshold
#  correction_method = "g_SCS", #default multiple testing correction method to reduce false positives
#  domain_scope = "annotated", #only consider annoted genes
#  custom_bg = NULL, 
#  numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE
#)
#head(gost_res$result)
#gostplot(gost_res2, capped = TRUE, interactive = FALSE)
#publish_gosttable(gost_res)


#enrichment analysis on differential expression ms vs healthy control
#show all
gost_res2 <- gost(
  query = unlist(diff_expr_df[diff_expr_df$padj<0.05, 'Gene']), 
  significant=FALSE, 
  organism = "hsapiens"
  ) 
gostplot(gost_res2, capped=TRUE, interactive=FALSE)
head(gost_res2$results)
#significant only
gost_res3 <- gost(
  query = unlist(diff_expr_df[diff_expr_df$padj<0.05, 'Gene']), 
  significant=TRUE, 
  organism = "hsapiens"
  ) 
png("SRP164913_gprofiler_gostplot.png")
gostplot(gost_res3, capped=TRUE, interactive=FALSE)
dev.off()
head(gost_res3$results)

#head of both returns
publish_gosttable(gost_res3)

