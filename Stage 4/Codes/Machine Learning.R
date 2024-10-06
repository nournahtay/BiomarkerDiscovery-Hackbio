#Load bioconductor and TCGAbiolinks
library("BiocManager")
library("TCGAbiolinks")
library("dplyr")
library("limma")
library("edgeR")
library("caret")
library("gplots")
library("sesameData")
library("SummarizedExperiment")

#obtain tcga_lgg data
 
 LGG_g <- GDCquery(project = "TCGA-LGG",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification")
 GDCdownload(LGG_g)
 LGG_exp.Data <- GDCprepare(LGG_g)
 LGG.metaData <- data.frame("IDH_status" = LGG_exp.Data$paper_IDH.status,
                            "Barcode" = LGG_exp.Data$barcode)
 
#df to explore the data and subset important metrics
LGG_df <- data.frame("Mutation_count" = LGG_exp.Data@colData@listData[["paper_Mutation.Count"]],
                     "Barcode" = LGG_exp.Data$barcode,
                     "Histology" = LGG_exp.Data@colData@listData[["paper_Histology"]],
                     "tumor_type" = LGG_exp.Data@colData@listData[["tumor_descriptor"]],
                     "Gender" = LGG_exp.Data@colData@listData[["gender"]],
                     "Grade" = LGG_exp.Data@colData@listData[["paper_Grade"]],
                     "MGMTpromoterStatus" = LGG_exp.Data@colData@listData[["paper_MGMT.promoter.status"]],
                     "Expressioncluster" = LGG_exp.Data@colData@listData[["paper_Pan.Glioma.RNA.Expression.Cluster"]])

#Stranded data was selected                      
LGGRaw <- assays(LGG_exp.Data)
dim(LGGRaw$stranded_first)
SelectedBarcodes<- LGG.metaData$Barcode
                                        
SelectedData <- LGGRaw$stranded_first[, c(SelectedBarcodes)]
#Variance thresholding
high_variance_cols <- nearZeroVar(SelectedData, saveMetrics = TRUE)
#filter
lgg_filtered <- SelectedData[, !high_variance_cols$nzv]
#normalise
lgg_normalized <- TCGAanalyze_Normalization(tabDF = lgg_filtered, geneInfo = geneInfoHT, method = "geneLength")

#randomn walk for feature selection
 
install.packages("igraph")
library(igraph)
lgg_gene_correlation <- cor(lgg_normalized)
# Compute correlations
lgg_similarity_graph <- graph_from_adjacency_matrix(lgg_gene_correlation, mode = "undirected", weighted = TRUE)
set.seed(123)  # Set seed for reproducibility
start_vertex <- sample(V(lgg_similarity_graph), 1)

lgg_walk <- random_walk(lgg_similarity_graph, start = start_vertex, steps = 100)
gene_ranking <- order(lgg_walk, decreasing=TRUE)  # Rank genes based on importance

top_n <- 100
selected_genes <- gene_ranking[1:top_n]  # Choose top 'n' genes
library(kknn)
library(randomForest)
selected_features <- lgg_normalized[, selected_genes] # Subset selected genes from Random Walk
View(selected_features)
View(lgg_normalized)

selected_features.t <- t(selected_features)  # Transpose the training data

library(class)
library(cluster)


kmeans_result <- kmeans(selected_features.t, centers = 4 )  # Unsupervised clustering
expressioncluster <- kmeans_result$cluster
expressioncluster <- as.factor(expressioncluster)


#visualize the clusters
library(ggplot2)
cluster_df <- data.frame(PCA1 = kmeans_result$centers[, 1], PCA2 = kmeans_result$centers[, 2], Cluster = factor(1:4))
ggplot(cluster_df, aes(x = PCA1, y = PCA2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering Results")

#Randomn forest
library(randomForest)
set.seed(20)  # Set seed for reproducibility
rf_model <- randomForest(selected_features.t, y = expressioncluster, ntree = 100)  # Train the model
print(rf_model)

#assess feature importance

# Calculate feature importance
importance_values <- importance(rf_model)

# Display the importance values
print(importance_values)

# Optionally, visualize the importance
varImpPlot(rf_model)

