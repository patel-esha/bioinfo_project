#Hannah Luft Assignment 3

#file paths
data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
diff_expr_rslts_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")
#load in the expression data and metadata
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Symbol")
metadata <- readr::read_tsv(metadata_file)
#load in differential expression data
diff_expr <- readr::read_tsv(diff_expr_rslts_file)

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


#libraries
library(tidyclust)
library(tidyverse)
#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(DESeq2)
library(ggplot2)
library(magrittr)
library(ComplexHeatmap)
#install.packages("ggalluvial")
library(ggalluvial)
#install.packages("flextable")
set.seed(1234)


#subset data into top 5000 differentially expressed genes
#remove NA padj values (correspond to pvals of 1)
#diff_expr_nona <- na.omit(diff_expr)
#diff_expr_ordered <- diff_expr_nona[order(diff_expr_nona$padj),] #sort so smallest padj to largest
#only has 1980 rows after removing all NAs. 
#dim(diff_expr_ordered) 
#list of gene names in order that are the top differentially expressed
#top_genes <- diff_expr[1:1980,1] 
#get only those genes from the expression data
#diff_genes_expr_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes$Gene,]
#make the row names the gene symbols
#diff_genes_expr_df <- diff_genes_expr_df %>% remove_rownames %>% column_to_rownames(var="Symbol")

# log scale
#epsilon <- 1e-6
#culled_expression_df <- log2(culled_expression_df + epsilon)
# order genes by variance
culled_expression_df$variance <- apply(culled_expression_df, 1, var) 
exp_ordered <- culled_expression_df[order(culled_expression_df$variance, decreasing = TRUE), ]
expressions <- select(exp_ordered, -variance)
exp_top10 <- expressions[1:10, ]
exp_top100 <- expressions[1:100, ]
exp_top1000 <- expressions[1:1000, ]
exp_top5000 <- expressions[1:5000, ]
exp_top10000 <- expressions[1:10000, ]
diff_genes_expr_df <- exp_top5000

#Concensus cluster plus algorithm
d = data.matrix(diff_genes_expr_df)
#run concensus cluster plus
results_path <- file.path(results_dir, "consensusClusterPlusResults/k=6")
results <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=results_path,distance="pearson",plot="png")
#cluster and item consensus
icl = calcICL(results,title=results_path,plot="png")


#alter k
#k=10
results_path <- file.path(results_dir, "consensusClusterPlusResults/k=10")
results2 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,distance="pearson",plot="png")
icl2 = calcICL(results2,title=results_path,plot="png")

#k=20
results_path <- file.path(results_dir, "consensusClusterPlusResults/k=20")
results3 <- ConsensusClusterPlus(d,maxK=20,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl3 = calcICL(results3,title=results_path,plot="png")

#rerun with dif num of genes
#10 genes
#top_genes_10 <- diff_expr[1:10,1]
#genes_expr_10_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_10$Gene,]
#genes_expr_10_df <- genes_expr_10_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d4 = data.matrix(exp_top10)
results_path <- file.path(results_dir, "consensusClusterPlusResults/g=10")
results4 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl4 = calcICL(results4,title=results_path,plot="png")

#100 genes
#top_genes_100 <- diff_expr[1:100,1]
#genes_expr_100_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_100$Gene,]
#genes_expr_10_df <- genes_expr_100_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d5 = data.matrix(exp_top100)
results_path <- file.path(results_dir, "consensusClusterPlusResults/g=100")
results5 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl5 = calcICL(results5,title=results_path,plot="png")

#1000 genes
#top_genes_1000 <- diff_expr[1:1000,1]
#genes_expr_1000_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_1000$Gene,]
#genes_expr_1000_df <- genes_expr_1000_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d6 = data.matrix(exp_top1000)
results_path <- file.path(results_dir, "consensusClusterPlusResults/g=1000")
results6 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl6 = calcICL(results6,title=results_path,plot="png")

#chi squared test with each pair of clustering results
#turn data into a contingency table(s): 10vs100, 10vs1000, 10vs1980, 100vs1000, 100vs1980, 1000vs1980
#           g=10    g=100
# cluster1
# cluster2
#prepare data
clusterAssignmentsk2g1980 <- results[[2]]["consensusClass"]
clusterAssignmentsk2g10 <- results4[[2]]["consensusClass"]
clusterAssignmentsk2g100 <- results5[[2]]["consensusClass"]
clusterAssignmentsk2g1000 <- results6[[2]]["consensusClass"]

chi10vs100 <- data.frame(g10 = c(unname(table(clusterAssignmentsk2g10[[1]])[1]),unname(table(clusterAssignmentsk2g10[[1]])[2])),
                         g100 = c(unname(table(clusterAssignmentsk2g100[[1]])[1]), unname(table(clusterAssignmentsk2g100[[1]])[2])))
chi10vs1000 <- data.frame(g10 = c(unname(table(clusterAssignmentsk2g10[[1]])[1]),unname(table(clusterAssignmentsk2g10[[1]])[2])),
                         g1000 = c(unname(table(clusterAssignmentsk2g1000[[1]])[1]), unname(table(clusterAssignmentsk2g1000[[1]])[2])))
chi10vs1980 <- data.frame(g10 = c(unname(table(clusterAssignmentsk2g10[[1]])[1]),unname(table(clusterAssignmentsk2g10[[1]])[2])),
                         g1980 = c(unname(table(clusterAssignmentsk2g1980[[1]])[1]), unname(table(clusterAssignmentsk2g1980[[1]])[2])))
chi100vs1000 <- data.frame(g100 = c(unname(table(clusterAssignmentsk2g100[[1]])[1]),unname(table(clusterAssignmentsk2g100[[1]])[2])),
                         g1000 = c(unname(table(clusterAssignmentsk2g1000[[1]])[1]), unname(table(clusterAssignmentsk2g1000[[1]])[2])))
chi100vs1980 <- data.frame(g100 = c(unname(table(clusterAssignmentsk2g100[[1]])[1]),unname(table(clusterAssignmentsk2g100[[1]])[2])),
                         g1980 = c(unname(table(clusterAssignmentsk2g1980[[1]])[1]), unname(table(clusterAssignmentsk2g1980[[1]])[2])))
chi1000vs1980 <- data.frame(g1000 = c(unname(table(clusterAssignmentsk2g1000[[1]])[1]),unname(table(clusterAssignmentsk2g1000[[1]])[2])),
                         g1980 = c(unname(table(clusterAssignmentsk2g1980[[1]])[1]), unname(table(clusterAssignmentsk2g1980[[1]])[2])))

#run the chi squared test on tables
chi_10v100_2k <- chisq.test(chi10vs100)
chi_10v1000_2k <- chisq.test(chi10vs1000)
chi_10v1980_2k <- chisq.test(chi10vs1980)
chi_100v1000_2k <- chisq.test(chi100vs1000)
chi_100v1980_2k <- chisq.test(chi100vs1980)
chi_1000v1980_2k <- chisq.test(chi1000vs1980)
#all have x-squared value of 0 and p value of 1

#do it again with k=3
clusterAssignmentsk3g1980 <- results[[3]]["consensusClass"]
clusterAssignmentsk3g10 <- results4[[3]]["consensusClass"]
clusterAssignmentsk3g100 <- results5[[3]]["consensusClass"]
clusterAssignmentsk3g1000 <- results6[[3]]["consensusClass"]
#           g10  g100
#cluster 1
#cluster 2
#cluster 3

chi3k10vs100 <- data.frame(g10 = c(unname(table(clusterAssignmentsk3g10[[1]])[1]),unname(table(clusterAssignmentsk3g10[[1]])[2]), unname(table(clusterAssignmentsk3g10[[1]])[3])),
                         g100 = c(unname(table(clusterAssignmentsk3g100[[1]])[1]), unname(table(clusterAssignmentsk3g100[[1]])[2]), unname(table(clusterAssignmentsk3g100[[1]])[3])))
chi3k10vs1000 <- data.frame(g10 = c(unname(table(clusterAssignmentsk3g10[[1]])[1]),unname(table(clusterAssignmentsk3g10[[1]])[2]), unname(table(clusterAssignmentsk3g10[[1]])[3])),
                           g1000 = c(unname(table(clusterAssignmentsk3g1000[[1]])[1]), unname(table(clusterAssignmentsk3g1000[[1]])[2]), unname(table(clusterAssignmentsk3g1000[[1]])[3])))
chi3k10vs1980 <- data.frame(g10 = c(unname(table(clusterAssignmentsk3g10[[1]])[1]),unname(table(clusterAssignmentsk3g10[[1]])[2]), unname(table(clusterAssignmentsk3g10[[1]])[3])),
                           g1980 = c(unname(table(clusterAssignmentsk3g1980[[1]])[1]), unname(table(clusterAssignmentsk3g1980[[1]])[2]), unname(table(clusterAssignmentsk3g1980[[1]])[3])))
chi3k100vs1000 <- data.frame(g100 = c(unname(table(clusterAssignmentsk3g100[[1]])[1]),unname(table(clusterAssignmentsk3g100[[1]])[2]), unname(table(clusterAssignmentsk3g100[[1]])[3])),
                           g1000 = c(unname(table(clusterAssignmentsk3g1000[[1]])[1]), unname(table(clusterAssignmentsk3g1000[[1]])[2]), unname(table(clusterAssignmentsk3g1000[[1]])[3])))
chi3k100vs1980 <- data.frame(g100 = c(unname(table(clusterAssignmentsk3g100[[1]])[1]),unname(table(clusterAssignmentsk3g100[[1]])[2]), unname(table(clusterAssignmentsk3g100[[1]])[3])),
                           g1980 = c(unname(table(clusterAssignmentsk3g1980[[1]])[1]), unname(table(clusterAssignmentsk3g1980[[1]])[2]), unname(table(clusterAssignmentsk3g1980[[1]])[3])))
chi3k1000vs1980 <- data.frame(g1000 = c(unname(table(clusterAssignmentsk3g1000[[1]])[1]),unname(table(clusterAssignmentsk3g1000[[1]])[2]), unname(table(clusterAssignmentsk3g1000[[1]])[3])),
                           g1980 = c(unname(table(clusterAssignmentsk3g1980[[1]])[1]), unname(table(clusterAssignmentsk3g1980[[1]])[2]), unname(table(clusterAssignmentsk3g1980[[1]])[3])))
chi_10v100_3k <- chisq.test(chi3k10vs100)
chi_10v1000_3k <- chisq.test(chi3k10vs1000)
chi_10v1980_3k <- chisq.test(chi3k10vs1980)
chi_100v1000_3k <- chisq.test(chi3k100vs1000)
chi_100v1980_3k <- chisq.test(chi3k100vs1980)
chi_1000v1980_3k <- chisq.test(chi3k1000vs1980)
#all have x-squared value of 0 and p value of 1

#plot to see what is going on- k=2 clusters
metadata_clustered_2k <- metadata
metadata_clustered_2k$cluster_2k <- clusterAssignmentsk2g1980[[1]]
metadata_clustered_2k <- metadata_clustered_2k %>%
  dplyr::mutate(
    cluster_2k = factor(
      cluster_2k,
      levels = c("1", "2")
    )
  )
filtered_data_2k_df <- exp_top5000 %>% dplyr::filter(rowSums(.) >=1)
filtered_data_2k_df <- round(filtered_data_2k_df)
dds_2k <- DESeqDataSetFromMatrix(
  countData = filtered_data_2k_df,
  colData = metadata_clustered_2k,
  design = ~1
)
orows <- sum( rowMeans( counts(dds_2k, normalized=FALSE)) > 5 )
dds_norm_2k <- vst(dds_2k, nsub=20)
plotPCA(
  dds_norm_2k,
  intgroup = "cluster_2k"
)
pca_results_2k <-
  plotPCA(
    dds_norm_2k,
    intgroup = c("cluster_2k"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
annotated_pca_plot_2k <- ggplot(
  pca_results_2k,
  aes(
    x = PC1,
    y = PC2,
    color = cluster_2k,
  )
) +
  geom_point()
# display plot
annotated_pca_plot_2k
results_path <- file.path(results_dir, "consensusClusterPlusResults")
ggsave(
  file.path(results_path, "consensusClusterPlusPlot2k.png"),
  plot = annotated_pca_plot_2k # the plot object that we want saved to file
)

#do it again with k=3 groups
metadata_clustered_3k <- metadata
metadata_clustered_3k$cluster_3k <- clusterAssignmentsk3g1980[[1]]
metadata_clustered_3k <- metadata_clustered_3k %>%
  dplyr::mutate(
    cluster_3k = factor(
      cluster_3k,
      levels = c("1", "2","3")
    )
  )
filtered_data_3k_df <- exp_top5000 %>% dplyr::filter(rowSums(.) >=1)
filtered_data_3k_df <- round(filtered_data_3k_df)
dds_3k <- DESeqDataSetFromMatrix(
  countData = filtered_data_3k_df,
  colData = metadata_clustered_3k,
  design = ~1
)
orows <- sum( rowMeans( counts(dds_3k, normalized=FALSE)) > 5 )
dds_norm_3k <- vst(dds_3k, nsub=20)
plotPCA(
  dds_norm_3k,
  intgroup = "cluster_3k"
)
pca_results_3k <-
  plotPCA(
    dds_norm_3k,
    intgroup = c("cluster_3k"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
annotated_pca_plot_3k <- ggplot(
  pca_results_3k,
  aes(
    x = PC1,
    y = PC2,
    color = cluster_3k,
  )
) +
  geom_point()
# display annotated plot
annotated_pca_plot_3k
results_path <- file.path(results_dir, "consensusClusterPlusResults")
ggsave(
  file.path(results_path, "consensusClusterPlusPlot3k.png"),
  plot = annotated_pca_plot_3k # the plot object that we want saved to file
)

#again with the 2 groups ms and hc we were trying to compare
metadata_clustered_groups <- metadata
metadata_clustered_groups <- metadata_clustered_groups %>%
  dplyr::mutate(
    refinebio_disease = factor(
      refinebio_disease,
      levels = c("hs", "ms")
    )
  )
pca_results_groups <-
  plotPCA(
    dds_norm_2k,
    intgroup = c("refinebio_disease"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
annotated_pca_plot_groups <- ggplot(
  pca_results_groups,
  aes(
    x = PC1,
    y = PC2,
    color = refinebio_disease,
  )
) +
  geom_point()
# display annotated plot
annotated_pca_plot_groups
results_path <- file.path(results_dir, "consensusClusterPlusResults")
ggsave(
  file.path(results_path, "consensusClusterPlusPlotgroups.png"),
  plot = annotated_pca_plot_groups # the plot object that we want saved to file
)


#alluvial diagram
#columns: cluster result (k=2-3)
#rows: samples
#create the dataframe
combined_metadata <- metadata
combined_metadata$cluster_2k <- clusterAssignmentsk2g1980[[1]]
combined_metadata$cluster_3k <-clusterAssignmentsk3g1980[[1]]
#group   k=2group  k=3group
#hc in group 1 of k=2 and group 1 of k=3
freq1 <- nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "1" &
                      cluster_3k == "1"
))
#hc in group 1 of k=2 and group 2 of k=3
freq2 <- nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "1" &
                      cluster_3k == "2"
))

#hc in group 1 of k=2 and group 3 of k=3
freq3 <-nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "1" &
                      cluster_3k == "3"
))
#hc in group 2 of k=2 and group 1 of k=3
freq4 <-nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "2" &
                      cluster_3k == "1"
))
#hc in group 2 of k=2 and group 2 of k=3
freq5 <-nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "2" &
                      cluster_3k == "2"
))
#hc in group 2 of k=2 and group 3 of k=3
freq6 <-nrow(subset(combined_metadata,
                    refinebio_disease == "hc" & 
                      cluster_2k == "2" &
                      cluster_3k == "3"
))

#ms in group 1 of k=2 and group 1 of k=3
freq7 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "1" &
                      cluster_3k == "1"
))
#ms in group 1 of k=2 and group 2 of k=3
freq8 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "1" &
                      cluster_3k == "2"
))
#ms in group 1 of k=2 and group 3 of k=3
freq9 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "1" &
                      cluster_3k == "3"
))
#ms in group 2 of k=2 and group 1 of k=3
freq10 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "2" &
                      cluster_3k == "1"
))
#ms in group 2 of k=2 and group 2 of k=3
freq11 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "2" &
                      cluster_3k == "2"
))
#ms in group 2 of k=2 and group 3 of k=3
freq12 <-nrow(subset(combined_metadata,
                    refinebio_disease == "ms" & 
                      cluster_2k == "2" &
                      cluster_3k == "3"
))

alluvial_df <- data.frame(
                disease_group =    c("hs", "hs", "hs", "hs", "hs", "hs", "ms", "ms", "ms", "ms", "ms", "ms"),
                cluster_k2_group = c("1",  "1",  "1",  "2",  "2",  "2",  "1",  "1",  "1",  "2",  "2",  "2"),
                cluster_k3_group = c("1",  "2",  "3",  "1",  "2",  "3",  "1",  "2",  "3",  "1",  "2",  "3"),
                freq = c(freq1, freq2, freq3, freq4, freq5, freq6, freq7, freq8, freq9, freq10, freq11, freq12)
              )

colors <- hcl.colors(2, "Red-Blue")
alluvial_graph <- ggplot(alluvial_df,
       aes(y = freq, axis1 = disease_group, axis2 = cluster_k2_group, axis3 = cluster_k3_group)) +
  geom_alluvium(aes(fill = disease_group), width = 1/12) +
  geom_stratum(width = 1/8) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Disease Group", "Cluster k=2", "Cluster k=3")) +
  scale_fill_manual(values = colors) +
  ggtitle("Cluster Groups by Healthy Control vs MS")

alluvial_graph
ggsave(
  file.path(results_path, "consensusClusterPlusAlluvialMap.png"),
  plot = alluvial_graph 
)

#heatmaps and dendrograms
#heatmap of the 1980 genes used. Add annotation side bars of clusters from each method and the two groups
#original 1980 genes: diff_genes_expr_df
#combined_metadata for the two groups, clusters k=2, and clusters k=3
column_ha = HeatmapAnnotation(groups = as.factor(combined_metadata$refinebio_disease), clusters2k = as.factor(combined_metadata$cluster_2k), clusters3k = as.factor(combined_metadata$cluster_3k))
#make the heatmap
png("results/consensusClusterPlusResults/consensusClusterPlusHeatmap.png")
ht <- Heatmap(
  as.matrix(exp_top1000), 
  name="Cluster Heatmap", 
  top_annotation = column_ha,
  column_title = "Samples", column_title_side = "bottom", 
  row_title = "Genes", row_title_side = "right" ,
  show_column_names = FALSE,
  show_row_names = FALSE
)
draw(ht)
dev.off()

#chi squared test of independence on two groups, for each clustering result
#chi table format: 
#         cluster1  cluster2
# group1
# group2 

#cluster k=2
chi_disease_2k <- data.frame(
  cluster1 = c(
                nrow(subset(combined_metadata, refinebio_disease == "hc" & cluster_2k == "1")),
                nrow(subset(combined_metadata, refinebio_disease == "ms" & cluster_2k == "1"))
              ),
  cluster2 = c(
                nrow(subset(combined_metadata, refinebio_disease == "hc" & cluster_2k == "2")),
                nrow(subset(combined_metadata, refinebio_disease == "ms" & cluster_2k == "2"))
  )
)
chi_groupsv2k <- chisq.test(chi_disease_2k)

#cluster k=3
#         cluster1  cluster2  cluster3
# group1
# group2 
chi_disease_3k <- data.frame(
  cluster1 = c(
    nrow(subset(combined_metadata, refinebio_disease == "hc" & cluster_3k == "1")),
    nrow(subset(combined_metadata, refinebio_disease == "ms" & cluster_3k == "1"))
  ),
  cluster2 = c(
    nrow(subset(combined_metadata, refinebio_disease == "hc" & cluster_3k == "2")),
    nrow(subset(combined_metadata, refinebio_disease == "ms" & cluster_3k == "2"))
  ),
  cluster3 = c(
    nrow(subset(combined_metadata, refinebio_disease == "hc" & cluster_3k == "3")),
    nrow(subset(combined_metadata, refinebio_disease == "ms" & cluster_3k == "3"))
  )
)
chi_groupsv3k <- chisq.test(chi_disease_3k)


#compile all chi square test results into one table
chi_table_df <- data.frame(
  Test = c("chi_10v100_2k", "chi_10v1000_2k", "chi_10v1980_2k", "chi_100v1000_2k", "chi_100v1980_2k", "chi_1000v1980_2k", 
           "chi_10v100_3k", "chi_10v1000_3k", "chi_10v1980_3k", "chi_100v1000_3k", "chi_100v1980_3k", "chi_1000v1980_3k",
           "chi_groupsv2k", "chi_groupsv3k"),
  X_Squared_Value = c(unname(chi_10v100_2k$statistic), unname(chi_10v1000_2k$statistic), unname(chi_10v1980_2k$statistic), unname(chi_100v1000_2k$statistic), unname(chi_100v1980_2k$statistic), unname(chi_1000v1980_2k$statistic), 
                      unname(chi_10v100_3k$statistic), unname(chi_10v1000_3k$statistic), unname(chi_10v1980_3k$statistic), unname(chi_100v1000_3k$statistic), unname(chi_100v1980_3k$statistic), unname(chi_1000v1980_3k$statistic),
                      unname(chi_groupsv2k$statistic), unname(chi_groupsv3k$statistic)),
  p_Value = c(chi_10v100_2k$p.value, chi_10v1000_2k$p.value, chi_10v1980_2k$p.value, chi_100v1000_2k$p.value, chi_100v1980_2k$p.value, chi_1000v1980_2k$p.value, 
              chi_10v100_3k$p.value, chi_10v1000_3k$p.value, chi_10v1980_3k$p.value, chi_100v1000_3k$p.value, chi_100v1980_3k$p.value, chi_1000v1980_3k$p.value,
              chi_groupsv2k$p.value, chi_groupsv3k$p.value)
)
#padjust all stat test results for multiple hypothesis testing
p_list <- as.list(chi_table_df$p_Value)
p_adj_list <- p.adjust(p_list, method="bonferroni")
chi_table_df$p_adj <- p_adj_list

#save the table to a file
fileName <- file.path(results_path, "ConsensusClusterPlus_X-Squared_Table.csv")
write_csv(chi_table_df, fileName)

