#Hannah Luft Assignment 3

data_dir <- file.path("data", "SRP164913")
data_file <- file.path(data_dir, "SRP164913_HUGO.tsv")
metadata_file <- file.path(data_dir, "metadata_SRP164913.tsv")
results_dir <- file.path("results")
plots_dir <- file.path("plots")
diff_expr_rslts_file <- file.path(results_dir, "SRP164913_diff_expr_results.tsv")
#load in the expression data and metadata
expression_df <- readr::read_tsv(data_file)
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
set.seed(1234)


#subset data into top 5000 differentially expressed genes
#remove NA padj values (correspond to pvals of 1)
diff_expr_nona <- na.omit(diff_expr)
diff_expr_ordered <- diff_expr_nona[order(diff_expr_nona$padj),] #sort so smallest padj to largest
dim(diff_expr_ordered) #only has 1980 rows after removing all NAs. 
top_genes <- diff_expr[1:1980,1] #list of gene names in order that are the top differentially expressed
#get only those genes from the expression data
diff_genes_expr_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes$Gene,]
#make the row names the gene symbols
diff_genes_expr_df <- diff_genes_expr_df %>% remove_rownames %>% column_to_rownames(var="Symbol")


#Concensus cluster plus algorithm
d = data.matrix(diff_genes_expr_df)
#run concensus cluster plus
results_path <- file.path(results_dir, "consensusClusterPlusResults/k=6")
results <- ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=results_path,distance="pearson",plot="png")
#cluster and item consensus
icl = calcICL(results,title=results_path,plot="png")


#alter k
results_path <- file.path(results_dir, "consensusClusterPlusResults/k=10")
results2 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,distance="pearson",plot="png")
icl2 = calcICL(results2,title=results_path,plot="png")

results_path <- file.path(results_dir, "consensusClusterPlusResults/k=20")
results3 <- ConsensusClusterPlus(d,maxK=20,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl3 = calcICL(results3,title=results_path,plot="png")

#rerun with dif num of genes
top_genes_10 <- diff_expr[1:10,1]
genes_expr_10_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_10$Gene,]
genes_expr_10_df <- genes_expr_10_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d4 = data.matrix(genes_expr_10_df)
results_path <- file.path(results_dir, "consensusClusterPlusResults/g=10")
results4 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl4 = calcICL(results4,title=results_path,plot="png")

top_genes_100 <- diff_expr[1:100,1]
genes_expr_100_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_100$Gene,]
genes_expr_10_df <- genes_expr_100_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d5 = data.matrix(genes_expr_100_df)
results_path <- file.path(results_dir, "consensusClusterPlusResults/g=100")
results5 <- ConsensusClusterPlus(d,maxK=10,reps=50,pItem=0.8,pFeature=1,title=results_path,clusterAlg="hc",distance="pearson",plot="png")
icl5 = calcICL(results5,title=results_path,plot="png")

top_genes_1000 <- diff_expr[1:1000,1]
genes_expr_1000_df <- culled_expression_df[culled_expression_df$Symbol %in% top_genes_1000$Gene,]
genes_expr_1000_df <- genes_expr_1000_df %>% remove_rownames %>% column_to_rownames(var="Symbol")
d6 = data.matrix(genes_expr_1000_df)
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
chisq.test(chi10vs100)
chisq.test(chi10vs1000)
chisq.test(chi10vs1980)
chisq.test(chi100vs1000)
chisq.test(chi100vs1980)
chisq.test(chi1000vs1980)
#all have x-squared value of 0 and p value of 1

#do it with k=3
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
chisq.test(chi3k10vs100)
chisq.test(chi3k10vs1000)
chisq.test(chi3k10vs1980)
chisq.test(chi3k100vs1000)
chisq.test(chi3k100vs1980)
chisq.test(chi3k1000vs1980)
#all have x-squared value of 0 and p value of 1

#plot
metadata_clustered <- metadata
metadata_clustered$cluster <- clusterAssignmentsk2g1980[[1]]
metadata_clustered <- metadata_clustered %>%
  dplyr::mutate(
    cluster = factor(
      cluster,
      # specify the possible levels in the order we want them to appear
      levels = c("1", "2")
    )
  )
filtered_data_df <- diff_genes_expr_df %>% dplyr::filter(rowSums(.) >=5)
filtered_data_df <- round(filtered_data_df)
dds <- DESeqDataSetFromMatrix(
  countData = filtered_data_df,
  colData = metadata_clustered,
  design = ~1
)
orows <- sum( rowMeans( counts(dds, normalized=FALSE)) > 5 )
dds_norm <- vst(dds, nsub=20)
plotPCA(
  dds_norm,
  intgroup = "cluster"
)
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("cluster"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_disease` group
    color = cluster,
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
annotated_pca_plot
results_path <- file.path(results_dir, "consensusClusterPlusResults")
ggsave(
  file.path(results_path, "consensusClusterPlusPlot2k.png"),
  plot = annotated_pca_plot # the plot object that we want saved to file
)

#do it again with 3 groups
metadata_clustered <- metadata
metadata_clustered$cluster <- clusterAssignmentsk3g1980[[1]]
metadata_clustered <- metadata_clustered %>%
  dplyr::mutate(
    cluster = factor(
      cluster,
      # specify the possible levels in the order we want them to appear
      levels = c("1", "2","3")
    )
  )
filtered_data_df <- diff_genes_expr_df %>% dplyr::filter(rowSums(.) >=5)
filtered_data_df <- round(filtered_data_df)
dds <- DESeqDataSetFromMatrix(
  countData = filtered_data_df,
  colData = metadata_clustered,
  design = ~1
)
orows <- sum( rowMeans( counts(dds, normalized=FALSE)) > 5 )
dds_norm <- vst(dds, nsub=20)
plotPCA(
  dds_norm,
  intgroup = "cluster"
)
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("cluster"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_disease` group
    color = cluster,
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
annotated_pca_plot
results_path <- file.path(results_dir, "consensusClusterPlusResults")
ggsave(
  file.path(results_path, "consensusClusterPlusPlot3k.png"),
  plot = annotated_pca_plot # the plot object that we want saved to file
)



#alluvial diagram


#heatmaps and dendrograms
#heatmap of the 1980 genes used. Add annotation side bars of clusters from each method and the two groups

#chi squared test of independence on two groups, for each clustering result
#chi table format: 
#         cluster1  cluster2
# group1
# group2 


#padjust all stat test results for multiple hypothesis testing




