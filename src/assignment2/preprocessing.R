"Data Preparation - Convert the Gene IDs from Ensemble to Hugo (Symbol).
The resulting file is written to 'data/SRP164913_HUGO.tsv'.

Make sure to run 'dependencies.rmd' before

This Script is based on the following tutorial:
https://alexslemonade.github.io/refinebio-examples/03-rnaseq/gene-id-annotation_rnaseq_01_ensembl.html
"

# Attach the library
library(org.Hs.eg.db)
library(AnnotationDbi) 

# We will need this so we can use the pipe: %>%
library(magrittr)

# Install ggplot2
library(ggplot2)


# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define the file path
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")

# check if files exists
file.exists(data_file)
file.exists(metadata_file)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated Symbol IDs
mapped_list <- mapIds(
  org.Hs.eg.db, 
  keys = expression_df$Gene,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list"
)
# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)

summary(as.factor(mapped_df$Symbol), maxsum = 10)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "symbol_id_count") %>%
  # Arrange by the genes with the highest number of Symbol IDs mapped
  dplyr::arrange(desc(symbol_id_count))

collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the Symbol IDs `mapped_df` into one column named `all_symbol_ids`
  dplyr::summarize(all_symbol_ids = paste(Symbol, collapse = ";"))

collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_symbol_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbol_ids, ";")) %>%
  # We only need a preview here
  head()

# keep only the first mapping
merge_mapped_df <- data.frame(
  "Symbol" = mapIds(
    org.Hs.eg.db,
    keys = expression_df$Gene,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Remove all rows where Symbol is empty
  dplyr::filter(!is.na(Symbol) & Symbol != "") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))

# Aggregate all duplicated symbols by taking the median of their expression values
final_merged_df <- aggregate(. ~ Symbol, data = merge_mapped_df, FUN = median, na.rm = TRUE)

# Set Symbol as identifier column
rownames(final_merged_df) <- final_merged_df$Symbol
head(final_merged_df)
final <- subset(final_merged_df, select = -c(Ensembl, all_symbol_ids))
head(final)
# Write mapped and annotated data frame to output file
readr::write_tsv(final, file.path(
  data_dir,
  "SRP164913_HUGO.tsv"
))

