# Define the file paths 
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")

# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)
# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)

#install packages
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

#attach libs
# Attach the DESeq2 library
library(DESeq2)
# Attach the ggplot2 library for plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)
set.seed(12345)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)
# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Symbol")
# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)
# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)
head(metadata$refinebio_title)

#extract metadata
metadata <- metadata %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(disease_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "HC") ~ "healthy",
    stringr::str_detect(refinebio_title, "MS") ~ "ms", 
    stringr::str_detect(refinebio_title, "HAM") ~ "ham/tsp"
  ))
# Let's take a look at the original metadata column's info
# and our new `mutation_status` column
dplyr::select(metadata, refinebio_title, disease_status)
# Print out a preview of `mutation_status`
str(metadata$disease_status)

# Make mutation_status a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    disease_status = factor(disease_status, levels = c("healthy", "ms", "ham/tsp"))
  )

levels(metadata$disease_status)

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)
# round all expression counts
gene_matrix <- round(filtered_expression_df)
ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~disease_status
)

deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 3, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)

# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)
plotCounts(ddset, gene = "LINC01317", intgroup = "disease_status")

readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "SRP164913_diff_expr_results.tsv" 
  )
)

#get top 50 to save as file
top_50 = top_n(deseq_df, 50, wt=log2FoldChange)
head(top_50)
#save to file
readr::write_tsv(
  top_50,
  file.path(
    results_dir, 
    "SRP164913_diff_expr_top_50_results.tsv"
  )
)


# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.001 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
volcano_plot
#volcano plot is VERY sparse, what's up with that? the p-vals are pretty high for most of them

ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "SRP164913_volcano_plot.png")
)
