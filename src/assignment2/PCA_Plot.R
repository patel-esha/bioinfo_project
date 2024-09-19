# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("data", "SRP164913")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
raw_data_file <- file.path(data_dir, "SRP164913.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")

results_dir <- file.path("results_dir")

hugo_data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")

plots_dir <- file.path("plots_dir")

# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(hugo_data_file)
# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)
# Read in data TSV file
expression_df <- readr::read_tsv(hugo_data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")


if (!("DESeq2" %in% installed.packages())) {
  # Install DESeq2
  BiocManager::install("DESeq2", update = FALSE)
}
# Attach the `DESeq2` library
library(DESeq2)
# Attach the `ggplot2` library for plotting
library(ggplot2)
# Set the seed so our results are reproducible:
set.seed(12345)
# Read in metadata TSV file (done earlier in code)
metadata <- readr::read_tsv(metadata_file)
# Read in data TSV file (done earlier)
expression_df <- readr::read_tsv(hugo_data_file) %>%
# Tuck away the gene ID column as row names, leaving only numeric values
  tibble::column_to_rownames("Gene")
# Make the sure the columns (samples) are in the same order as the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# convert the columns we will be using for annotation into factors
metadata <- metadata %>%
  dplyr::mutate(
    refinebio_disease = factor(
      refinebio_disease,
      # specify the possible levels in the order we want them to appear
      levels = c("ham/tsp", "hc", "ms")
    )
  )
#Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 5)
# The `DESeqDataSetFromMatrix()` function needs the values to be integers
filtered_expression_df <- round(filtered_expression_df)
# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df, # the counts values for all samples in our dataset
  colData = metadata, # annotation data for the samples in the counts data frame
  design = ~1 # Here we are not specifying a model
  # Replace with an appropriate design variable for your analysis
)
# Normalize and transform the data in the `DESeqDataSet` object
# using the `vst()` function from the `DESeq2` R package
dds_norm <- vst(dds)
plotPCA(
  dds_norm,
  intgroup = "refinebio_disease"
)
plotPCA(
  dds_norm,
  intgroup = c("refinebio_disease")
  # We are able to add another variable to the intgroup argument
  # by providing a vector of the variable names with `c()` function
)
# We first have to save the results of the `plotPCA()` function for use with `ggplot2`
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("refinebio_disease"),
    returnData = TRUE # This argument tells R to return the PCA values
  )

# Plot using `ggplot()` function and save to an object
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_treatment` group
    color = refinebio_disease,
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
annotated_pca_plot
# Save plot using `ggsave()` function
ggsave(
  file.path(plots_dir, "SRP164913_pca_plot.png"),
  # Replace with a file name relevant your plotted data
  plot = annotated_pca_plot # the plot object that we want saved to file
)


