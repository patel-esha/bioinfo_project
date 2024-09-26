#cluster profiler

data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
diff_expr_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")
plots_dir <- file.path("plots")


library(clusterProfiler)
library(readr)

expression_df <- readr::read_tsv(data_file)
gene_list <- as.list(expression_df[1])

diff_expr_df <- readr::read_tsv(diff_expr_file)

library(org.Hs.eg.db)  

# sig_genes contains the significantly differentially expressed genea
sig_genes <- diff_expr_df$Gene[diff_expr_df$padj < 0.05]

# Gene Ontology enrichment analysis
go_results <- enrichGO(gene = sig_genes,
                       OrgDb = org.Hs.eg.db,   
                       keyType = "SYMBOL",     
                       ont = "ALL",            
                       pvalueCutoff = 0.05,    
                       readable = TRUE)        

go_results_df <- as.data.frame(go_results)
output_file_tsv <- file.path(results_dir, "clusterprofiler_results.tsv")
write_tsv(go_results_df, file = output_file_tsv)
