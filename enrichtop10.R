#combining results
results_dir <- file.path("results")
topgo_results <- file.path(results_dir, "SRP164913_topGO_results.tsv")
clusterprofiler_results <- file.path(results_dir, "clusterprofiler_results.tsv")
wilcoxon_results <- file.path(results_dir, "Wilcoxon_GeneOnt/Enrichment.csv")

#read in files
topgo_df <- readr::read_tsv(topgo_results)
clusterprofiler_df <- readr::read_tsv(clusterprofiler_results)
wilcoxon_df <- readr::read_csv(wilcoxon_results)

head(topgo_df, 10)

clusterprofiler_sorted <- clusterprofiler_df[order(clusterprofiler_df$p.adjust,
                                decreasing = FALSE), ]
head(clusterprofiler_sorted, 10)

wilcoxon_sorted <- wilcoxon_df[order(clusterprofiler_df$p.adjust,
                                     decreasing = FALSE), ]
head(wilcoxon_df, 10)
