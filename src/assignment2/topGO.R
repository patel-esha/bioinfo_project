# Define the file paths 
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
diff_expr_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")
plots_dir <- file.path("plots")

# Load the required libraries
library("topGO")
all_genes <- read.delim(data_file, header = TRUE, stringsAsFactors = FALSE)
dGenes <- read.delim(diff_expr_file, header = TRUE, stringsAsFactors = FALSE)
head(dGenes)
# Create a named vector of p-values
pvals_list <- setNames(dGenes$pvalue, dGenes$Gene)

# Ensure all genes from the data file are included
geneList <- setNames(rep(1, length(all_genes$Gene)), all_genes$Gene)
geneList[names(pvals_list)] <- pvals_list

# Prepare the gene list for topGO by setting non-DE genes p-values to 1 (not significant)
geneList <- ifelse(is.na(geneList), 1, geneList)

# Create a function to select significant genes based on p-value threshold
pval_threshold <- 0.05
geneSelFun <- function(p) p < pval_threshold

# Now let's prepare the data for TopGO
tc <- new("topGOdata",
          ontology = "BP", # Change this if you want to use a different ontology
          allGenes = geneList,
          geneSel = geneSelFun,
          nodeSize = 10,
          annot = annFUN.org,
          mapping = "org.Hs.eg.db",
          ID = "symbol")

# Perform enrichment analysis using Fisher's exact test
resultFisher <- runTest(tc, algorithm = "classic", statistic = "fisher")

# To see most significant terms
allRes <- GenTable(tc, 
                   classicFisher = resultFisher, 
                   orderBy = "classicFisher", 
                   ranksOf = "classicFisher", 
                   topNodes = 10)
print(allRes)

