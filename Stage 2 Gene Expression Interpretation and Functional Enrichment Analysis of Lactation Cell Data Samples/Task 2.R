library(gplots)

url <- "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv"
gene_data <- read.csv(url, row.names = 1)

# show first few rows
head(gene_data)

# view data
# View(gene_data)
# generate heatmap
heatmap.2(as.matrix(gene_data), trace = 'none')

# scaling
heatmap.2(as.matrix(gene_data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE)

# adding colors (diverging colour palettes)
heatmap.2(as.matrix(gene_data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
            col=hcl.colors(100, palette = 'Blue-Red 3'))


dev.off()

# adding colors (sequential colour palettes)
heatmap.2(as.matrix(gene_data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'Purples 3'))

# getting column names
colnames(gene_data)

# selecting groups by index positions
group1 <- c(1,2,3,4,5)
group2 <- c(6,7,8,9,10)

# get groups 1 & 2 data (The groups behave similarly)
group1_data <- gene_data[, group1] #2005.01
group2_data <- gene_data[, group2] #1849.01

# get means
group1_mean <- rowMeans(group1_data)
group2_mean <- rowMeans(group2_data)

group1_mean
group2_mean

# get fold change
fold_change <- log2(group2_mean) - log2(group1_mean)  
fold_change

# get p-values
pvalues <- apply(gene_data, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
  
})

pvalues

# visualize the fold change and negative log of p-values
plot(fold_change, -log10(pvalues))

library(dplyr)

# Subset upregulated genes
upregulated_genes_1 <- group1_data %>%
  filter(fold_change > 2 & pvalues < 0.2)
upregulated_genes_2 <- group2_data %>%
  filter(fold_change > 2 & pvalues < 0.2)

# Subset downregulated genes
downregulated_genes_1 <- group1_data %>%
  filter(fold_change < -2 & pvalues < 0.2)
downregulated_genes_2 <- group2_data %>%
  filter(fold_change < -2 & pvalues < 0.2)

# Extract row names (gene names)
gene_namesup1 <- rownames(upregulated_genes_1)
gene_namesup2 <- rownames(upregulated_genes_2)
gene_namesdown1 <- rownames(downregulated_genes_1)
gene_namesdown2 <- rownames(downregulated_genes_2)

# Save the row names to a text file
write.table(gene_namesup1, file = "upregulated_genes_1.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(gene_namesup2, file = "upregulated_genes_2.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(gene_namesdown1, file = "downregulated_genes_1.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(gene_namesdown2, file = "downregulated_genes_2.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create data frame with pathways, number of genes, and FDR
top_pathways <- data.frame(
  Pathway = c("RNA Metabolism Regulation", 
              "RNA Biosynthesis Regulation", 
              "Regulation of transcription by RNA polymerase II", 
              "regulation of DNA-templated transcription", 
              "pattern specification process"),
  nGenes = c(13, 12, 12, 12, 6),  # Number of genes in each pathway
  FDR = c(0.0193, 0.0349, 0.0107, 0.0395, 0.0261)  # FDR values for each pathway
)

# Calculate -log10(FDR)
top_pathways$neg_log_fdr <- -log10(top_pathways$FDR)

library(ggplot2)

#Creating a Lollipop plor
ggplot(top_pathways, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_segment(aes(x = Pathway, xend = Pathway, y = 0, yend = nGenes), color = "grey") +
  geom_point(aes(size = neg_log_fdr), color = "blue") +
  labs(title = "Top 5 Pathways and Gene Count",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(FDR)") +
  theme_minimal() +
  coord_flip()

# Create a dot plot
ggplot(top_pathways, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_point(aes(size = neg_log_fdr), color = "blue") +
  labs(title = "Top 5 Pathways and Gene Count",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(FDR)") +
  theme_minimal() +
  coord_flip()

# Line plot with point size representing -log10(FDR)
ggplot(top_pathways, aes(x = reorder(Pathway, nGenes), y = nGenes, group = 1)) +
  geom_line(color = "darkgrey") +
  geom_point(aes(size = neg_log_fdr), color = "blue") +
  labs(title = "Top 5 Pathways and Gene Count",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(FDR)") +
  theme_minimal() +
  coord_flip()

# Bubble plot where point size and color represent -log10(FDR)
ggplot(top_pathways, aes(x = reorder(Pathway, nGenes), y = nGenes)) +
  geom_point(aes(size = neg_log_fdr, color = neg_log_fdr)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Top 5 Pathways and Gene Count",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(FDR)",
       color = "-log10(FDR)") +
  theme_minimal() +
  coord_flip()
