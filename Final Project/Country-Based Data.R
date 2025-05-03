# Load Important Packages
library("BiocManager")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("EDASeq")
library("gplots")
library("sesameData")
library("SummarizedExperiment")

# Get an overview on data
getProjectSummary("TCGA-STAD")

# Query the Data
STAD <- GDCquery(project = "TCGA-STAD",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

# Download the dataset
GDCdownload(STAD)

# Prepare the dataset
STAD.Data <- GDCprepare(STAD)

# Check the output of the data
head(STAD.Data)
View(STAD.Data)

# Explore the data
STAD.Data$gender
STAD.Data$paper_Country
STAD.Data$barcode
STAD.Data$patient
STAD.Data$paper_Country

#Check the Subsets of the metadata and add them in one data frame
STAD.META <- data.frame("Country" = STAD.Data$paper_Country,
                        "Stage" = STAD.Data$paper_TNM.Stage,
                        "Barcode" = STAD.Data$barcode)

#Omit Missing Values
STAD.META <- na.omit(STAD.META)
STAD.META <- subset(STAD.META, Stage != "X")

# Extract the raw data from the prepared dataset
STAD.Raw <- assays(STAD.Data)

#Select unstranded data
dim(STAD.Raw$unstranded)
View(STAD.Raw$unstranded)

# Group the stages into early and late stage
Early_stage <- c("Stage_IA", "Stage_IB", "Stage_IIA", "Stage_IIB")
Late_stage <- c("Stage_III", "Stage_IIIB", "Stage_IIIC", "Stage_IV")

# Subset the data
SelectedBarcodes <- c(subset(STAD.META, Country == "United_States" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "United_States" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Germany" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Germany" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Russia" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Russia" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Poland" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Poland" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Vietnam" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Vietnam" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Ukraine" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Ukraine" & Stage %in% Late_stage)$Barcode,
                      subset(STAD.META, Country == "Korea_South" & Stage %in% Early_stage)$Barcode,
                      subset(STAD.META, Country == "Korea_South" & Stage %in% Late_stage)$Barcode)

#Retrieve the unstranded data for the selected barcodes
SelectedData <- STAD.Raw$unstranded[, c(SelectedBarcodes)]
dim(SelectedData)
View(SelectedData)

#Normalize Data
normalized <- TCGAanalyze_Normalization(tabDF = SelectedData, geneInfo = geneInfoHT, method = "geneLength")

# Then Filter
filtered <- TCGAanalyze_Filtering(tabDF = normalized, method = "quantile", qnt.cut = 0.25)

# Create matrices for US samples
mat1.US.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "United_States" & Stage %in% Early_stage)$Barcode]
mat2.US.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "United_States" & Stage %in% Late_stage)$Barcode]

# Create matrices for Germany samples
mat1.Germany.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Germany" & Stage %in% Early_stage)$Barcode]
mat2.Germany.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Germany" & Stage %in% Late_stage)$Barcode]

# Create matrices for Russia samples
mat1.Russia.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Russia" & Stage %in% Early_stage)$Barcode]
mat2.Russia.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Russia" & Stage %in% Late_stage)$Barcode]

# Create matrices for Poland samples
mat1.Poland.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Poland" & Stage %in% Early_stage)$Barcode]
mat2.Poland.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Poland" & Stage %in% Late_stage)$Barcode]

# Create matrices for Vietnam samples
mat1.Vietnam.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Vietnam" & Stage %in% Early_stage)$Barcode]
mat2.Vietnam.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Vietnam" & Stage %in% Late_stage)$Barcode]

# Create matrices for Ukraine samples
mat1.Ukraine.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Ukraine" & Stage %in% Early_stage)$Barcode]
mat2.Ukraine.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Ukraine" & Stage %in% Late_stage)$Barcode]

# Create matrices for South Korea samples
mat1.Korea_South.Early <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Korea_South" & Stage %in% Early_stage)$Barcode]
mat2.Korea_South.Late <- filtered[, colnames(filtered) %in% subset(STAD.META, Country == "Korea_South" & Stage %in% Late_stage)$Barcode]

# Perform DEA for country Samples (early stages vs late stages)
results.US.1 <- TCGAanalyze_DEA(mat1 = mat1.US.Early, mat2 = mat2.US.Late,
                              Cond1type = "United_States %in% Early_stage",
                              Cond2type = "United_States %in% Late_stage",
                              pipeline = "edgeR",
                              fdr.cut = 0.1,
                              logFC.cut = 2)

results.Germany <- TCGAanalyze_DEA(mat1 = mat1.Germany.Early,
                              mat2 = mat2.Germany.Late,
                              Cond1type = "Germany %in% Early_stage",
                              Cond2type = "Germany %in% Late_stage",
                              pipeline = "edgeR",
                              fdr.cut = 0.1,
                              logFC.cut = 2)

results.Russia <- TCGAanalyze_DEA(mat1 = mat1.Russia.Early,
                                   mat2 = mat2.Russia.Late,
                                   Cond1type = "Russia %in% Early_stage",
                                   Cond2type = "Russia %in% Late_stage",
                                   pipeline = "edgeR",
                                  fdr.cut = 0.1,
                                  logFC.cut = 2)

results.Poland <- TCGAanalyze_DEA(mat1 = mat1.Poland.Early,
                                  mat2 = mat2.Poland.Late,
                                  Cond1type = "Poland %in% Early_stage",
                                  Cond2type = "Poland %in% Late_stage",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.1,
                                  logFC.cut = 2)

results.Vietnam <- TCGAanalyze_DEA(mat1 = mat1.Vietnam.Early,
                                  mat2 = mat2.Vietnam.Late,
                                  Cond1type = "Vietnam %in% Early_stage",
                                  Cond2type = "Vietnam %in% Late_stage",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.1,
                                  logFC.cut = 2)

results.Ukraine <- TCGAanalyze_DEA(mat1 = mat1.Ukraine.Early,
                                   mat2 = mat2.Ukraine.Late,
                                   Cond1type = "Ukraine %in% Early_stage",
                                   Cond2type = "Ukraine %in% Late_stage",
                                   pipeline = "edgeR",
                                   fdr.cut = 0.1,
                                   logFC.cut = 2)

results.Korea_South <- TCGAanalyze_DEA(mat1 = mat1.Korea_South.Early,
                                  mat2 = mat2.Korea_South.Late,
                                  Cond1type = "Korea_South %in% Early_stage",
                                  Cond2type = "Korea_South %in% Late_stage",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.1,
                                  logFC.cut = 2)

#Volcano Plots
US_volcano_plot <- plot(results.US$logFC, -log10( results.US$FDR))
Germany_volcano_plot <- plot(results.Germany$logFC, -log10(results.Germany$FDR))
Russia_volcano_plot <- plot(results.Russia$logFC, -log10(results.Russia$FDR))
Poland_volcano_plot <- plot(results.Poland$logFC, -log10(results.Poland$FDR))
Ukraine_volcano_plot <- plot(results.Ukraine$logFC, -log10(results.Ukraine$FDR))
Vietnam_volcano_plot <- plot(results.Vietnam$logFC, -log10(results.Vietnam$FDR))
Korea_South_volcano_plot <- plot(results.Korea_South$logFC, -log10(results.Korea_South$FDR))

#DEA with expression levels
US.level <- TCGAanalyze_LevelTab(results.US,
                                 "United_States %in% Early_stage",
                                 "United_States %in% Late_stage",
                                 mat1.US.Early, 
                                 mat2.US.Late)

Germany.level <- TCGAanalyze_LevelTab(results.Germany,
                                 "Germany %in% Early_stage",
                                 "Germany %in% Late_stage",
                                 mat1.Germany.Early, 
                                 mat2.Germany.Late)

Russia.level <- TCGAanalyze_LevelTab(results.Russia,
                                      "Russia %in% Early_stage",
                                      "Russia %in% Late_stage",
                                      mat1.Russia.Early, 
                                      mat2.Russia.Late)

Poland.level <- TCGAanalyze_LevelTab(results.Poland,
                                     "Poland %in% Early_stage",
                                     "Poland %in% Late_stage",
                                     mat1.Poland.Early, 
                                     mat2.Poland.Late)

Ukraine.level <- TCGAanalyze_LevelTab(results.Ukraine,
                                     "Ukraine %in% Early_stage",
                                     "Ukraine %in% Late_stage",
                                     mat1.Ukraine.Early, 
                                     mat2.Ukraine.Late)

Vietnam.level <- TCGAanalyze_LevelTab(results.Vietnam,
                                      "Vietnam %in% Early_stage",
                                      "Vietnam %in% Late_stage",
                                      mat1.Vietnam.Early, 
                                      mat2.Vietnam.Late)

Korea_South.level <- TCGAanalyze_LevelTab(results.Korea_South,
                                      "Korea_South %in% Early_stage",
                                      "Korea_South %in% Late_stage",
                                      mat1.Korea_South.Early, 
                                      mat2.Korea_South.Late)

#make sure to check your datasets

# Create vectors for the barcodes
US.Early.Barcodes <- subset(STAD.META, Country == "United_States" & Stage %in% Early_stage)$Barcode
US.Late.Barcodes <- subset(STAD.META, Country == "United_States" & Stage %in% Late_stage)$Barcode

Germany.Early.Barcodes <- subset(STAD.META, Country == "Germany" & Stage %in% Early_stage)$Barcode
Germany.Late.Barcodes <- subset(STAD.META, Country == "Germany" & Stage %in% Late_stage)$Barcode

Russia.Early.Barcodes <- subset(STAD.META, Country == "Russia" & Stage %in% Early_stage)$Barcode
Russia.Late.Barcodes <- subset(STAD.META, Country == "Russia" & Stage %in% Late_stage)$Barcode

Poland.Early.Barcodes <- subset(STAD.META, Country == "Poland" & Stage %in% Early_stage)$Barcode
Poland.Late.Barcodes <- subset(STAD.META, Country == "Poland" & Stage %in% Late_stage)$Barcode

Ukraine.Early.Barcodes <- subset(STAD.META, Country == "Ukraine" & Stage %in% Early_stage)$Barcode
Ukraine.Late.Barcodes <- subset(STAD.META, Country == "Ukraine" & Stage %in% Late_stage)$Barcode

Vietnam.Early.Barcodes <- subset(STAD.META, Country == "Vietnam" & Stage %in% Early_stage)$Barcode
Vietnam.Late.Barcodes <- subset(STAD.META, Country == "Vietnam" & Stage %in% Late_stage)$Barcode

Korea_South.Early.Barcodes <- subset(STAD.META, Country == "Korea_South" & Stage %in% Early_stage)$Barcode
Korea_South.Late.Barcodes <- subset(STAD.META, Country == "Korea_South" & Stage %in% Late_stage)$Barcode

# Combine the barcodes for heatmap data
Barcodes.US <- c(US.Early.Barcodes, US.Late.Barcodes)
Barcodes.Germany <- c(Germany.Early.Barcodes, Germany.Late.Barcodes)
Barcodes.Russia <- c(Russia.Early.Barcodes, Russia.Late.Barcodes)
Barcodes.Poland <- c(Poland.Early.Barcodes, Poland.Late.Barcodes)
Barcodes.Ukraine <- c(Ukraine.Early.Barcodes, Ukraine.Late.Barcodes)
Barcodes.Vietnam <- c(Vietnam.Early.Barcodes, Vietnam.Late.Barcodes)
Barcodes.Korea_South <- c(Korea_South.Early.Barcodes, Korea_South.Late.Barcodes)

# Subset the heatmap data 
heat.data.US <- filtered[rownames(US.level),
                              Barcodes.US]

heat.data.Germany <- filtered[rownames(Germany.level),
                         Barcodes.Germany]

heat.data.Russia <- filtered[rownames(Russia.level),
                         Barcodes.Russia]

heat.data.Poland <- filtered[rownames(Poland.level),
                         Barcodes.Poland]

heat.data.Ukraine <- filtered[rownames(Ukraine.level),
                         Barcodes.Ukraine]

heat.data.Vietnam <- filtered[rownames(Vietnam.level),
                         Barcodes.Vietnam]

heat.data.Korea_South <- filtered[rownames(Korea_South.level),
                         Barcodes.Korea_South]

# Create a color vector for the US samples
sample.US <- c(rep("midnightblue", 13), rep("red4", 6)) 

heatmap.2(x = as.matrix(heat.data.US),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of American Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.US)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Germany samples
sample.Germany <- c(rep("midnightblue", 11), rep("red4", 16)) 

heatmap.2(x = as.matrix(heat.data.Germany),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of German Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Germany)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Russia samples
sample.Russia <- c(rep("midnightblue", 39), rep("red4", 17)) 

heatmap.2(x = as.matrix(heat.data.Russia),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Russian Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Russia)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Poland samples
sample.Poland <- c(rep("midnightblue", 16), rep("red4", 10)) 

heatmap.2(x = as.matrix(heat.data.Poland),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Polish Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Poland)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Ukraine samples
sample.Ukraine <- c(rep("midnightblue", 16), rep("red4", 16)) 

heatmap.2(x = as.matrix(heat.data.Ukraine),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Ukrainian Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Ukraine)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the Vietnam samples
sample.Vietnam <- c(rep("midnightblue", 28), rep("red4", 8)) 

heatmap.2(x = as.matrix(heat.data.Vietnam),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Vietnamese Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Vietnam)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Create a color vector for the South Korean samples
sample.Korea_South <- c(rep("midnightblue", 18), rep("red4", 8)) 

heatmap.2(x = as.matrix(heat.data.Korea_South),
          col = hcl.colors(10, palette = 'Blue-Red 3'),  
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of South Korean Samples (Early Stage vs Late Stage)",
          na.color = 'black',
          ColSideColors = sample.Korea_South)

#Add a Legend
legend("topright", legend = c("Early stage", "Late stage"),
       fill = c("midnightblue", "red4"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

dev.off()

# Upregulated and downregulated genes for US samples
US.upreg.early <- rownames(subset(US.level, logFC > 1.5 & "United_States %in% Early_stage" > 0))
US.downreg.early <- rownames(subset(US.level, logFC < 1.5 & "United_States %in% Early_stage" > 0))
US.upreg.late <- rownames(subset(US.level, logFC > 1.5 & "United_States %in% Late_stage" > 0))
US.downreg.late <- rownames(subset(US.level, logFC < 1.5 & "United_States %in% Late_stage" > 0))


# Upregulated and downregulated genes for Germany samples
Germany.upreg.early <- rownames(subset(Germany.level, logFC > 1.5 & "Germany %in% Early_stage" > 0))
Germany.downreg.early <- rownames(subset(Germany.level, logFC < 1.5 & "Germany %in% Early_stage" > 0))
Germany.upreg.late <- rownames(subset(Germany.level, logFC > 1.5 & "Germany %in% Late_stage" > 0))
Germany.downreg.late <- rownames(subset(Germany.level, logFC < 1.5 & "Germany %in% Late_stage" > 0))

# Upregulated and downregulated genes for Russia samples
Russia.upreg.early <- rownames(subset(Russia.level, logFC > 1.5 & "Russia %in% Early_stage" > 0))
Russia.downreg.early <- rownames(subset(Russia.level, logFC < 1.5 & "Russia %in% Early_stage" > 0))
Russia.upreg.late <- rownames(subset(Russia.level, logFC > 1.5 & "Russia %in% Late_stage" > 0))
Russia.downreg.late <- rownames(subset(Russia.level, logFC < 1.5 & "Russia %in% Late_stage" > 0))

# Upregulated and downregulated genes for Poland samples
Poland.upreg.early <- rownames(subset(Poland.level, logFC > 1.5 & "Poland %in% Early_stage" > 0))
Poland.downreg.early <- rownames(subset(Poland.level, logFC < 1.5 & "Poland %in% Early_stage" > 0))
Poland.upreg.late <- rownames(subset(Poland.level, logFC > 1.5 & "Poland %in% Late_stage" > 0))
Poland.downreg.late <- rownames(subset(Poland.level, logFC < 1.5 & "Poland %in% Late_stage" > 0))

# Upregulated and downregulated genes for Ukraine samples
Ukraine.upreg.early <- rownames(subset(Ukraine.level, logFC > 2 & "Ukraine %in% Early_stage" > 0))
Ukraine.downreg.early <- rownames(subset(Ukraine.level, logFC < 2 & "Ukraine %in% Early_stage" > 0))
Ukraine.upreg.late <- rownames(subset(Ukraine.level, logFC > 2 & "Ukraine %in% Late_stage" > 0))
Ukraine.downreg.late <- rownames(subset(Ukraine.level, logFC < 2 & "Ukraine %in% Late_stage" > 0))

# Upregulated and downregulated genes for Vietnam samples
Vietnam.upreg.early <- rownames(subset(Vietnam.level, logFC > 2 & "Vietnam %in% Early_stage" > 0))
Vietnam.downreg.early <- rownames(subset(Vietnam.level, logFC < 2 & "Vietnam %in% Early_stage" > 0))
Vietnam.upreg.late <- rownames(subset(Vietnam.level, logFC > 2 & "Vietnam %in% Late_stage" > 0))
Vietnam.downreg.late <- rownames(subset(Vietnam.level, logFC < 2 & "Vietnam %in% Late_stage" > 0))

# Upregulated and downregulated genes for South Korea samples
Korea_South.upreg.early <- rownames(subset(Korea_South.level, logFC > 2 & "Korea_South %in% Early_stage" > 0))
Korea_South.downreg.early <- rownames(subset(Korea_South.level, logFC < 2 & "Korea_South %in% Early_stage" > 0))
Korea_South.upreg.late <- rownames(subset(Korea_South.level, logFC > 2 & "Korea_South %in% Late_stage" > 0))
Korea_South.downreg.late <- rownames(subset(Korea_South.level, logFC < 2 & "Korea_South %in% Late_stage" > 0))

library("biomaRt")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols for female samples
# Upregulated gene symbols for Early stage US samples
upregulated_Early_US_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = US.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage US samples
upregulated_Late_US_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = US.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage Germany samples
upregulated_Early_Germany_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Germany.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage Germany samples
upregulated_Late_Germany_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Germany.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage Russia samples
upregulated_Early_Russia_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Russia.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage Russia samples
upregulated_Late_Russia_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Russia.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage Poland samples
upregulated_Early_Poland_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Poland.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage Poland samples
upregulated_Late_Poland_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Poland.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage Ukraine samples
upregulated_Early_Ukraine_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Ukraine.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage Ukraine samples
upregulated_Late_Ukraine_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Ukraine.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage Vietnam samples
upregulated_Early_Vietnam_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Vietnam.upreg.early,
  mart = mart
)

# Upregulated gene symbols for Late stage Vietnam samples
upregulated_Late_Vietnam_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Vietnam.upreg.late,
  mart = mart
)

# Upregulated gene symbols for Early stage South Korea samples
upregulated_Early_Korea_South_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Korea_South.upreg.early,
  mart = mart
)


# Upregulated gene symbols for Late stage South Korea samples
upregulated_Late_Korea_South_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Korea_South.upreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage US samples
downregulated_Early_US_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = US.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage US samples
downregulated_Late_US_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = US.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage Germany samples
downregulated_Early_Germany_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Germany.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage Germany samples
downregulated_Late_Germany_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Germany.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage Russia samples
downregulated_Early_Russia_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Russia.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage Russia samples
downregulated_Late_Russia_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Russia.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage Poland samples
downregulated_Early_Poland_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Poland.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage Poland samples
downregulated_Late_Poland_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Poland.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage Ukraine samples
downregulated_Early_Ukraine_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Ukraine.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage Ukraine samples
downregulated_Late_Ukraine_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Ukraine.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage Vietnam samples
downregulated_Early_Vietnam_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Vietnam.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage Vietnam samples
downregulated_Late_Vietnam_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Vietnam.downreg.late,
  mart = mart
)

# Downregulated gene symbols for Early stage South Korea samples
downregulated_Early_Korea_South_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Korea_South.downreg.early,
  mart = mart
)

# Downregulated gene symbols for Late stage South Korea samples
downregulated_Late_Korea_South_symbols <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = Korea_South.downreg.late,
  mart = mart
)

# Enrichment Analysis
up.EA.US.early <- TCGAanalyze_EAcomplete(TFname = "United_States %in% Early_stage", RegulonList = upregulated_Early_US_symbols$hgnc_symbol)
up.EA.US.late <- TCGAanalyze_EAcomplete(TFname = "United_States %in% Late_stage", RegulonList = upregulated_Late_US_symbols$hgnc_symbol)
up.EA.Germany.early <- TCGAanalyze_EAcomplete(TFname = "Germany %in% Early_stage", RegulonList = upregulated_Early_Germany_symbols$hgnc_symbol)
up.EA.Germany.late <- TCGAanalyze_EAcomplete(TFname = "Germany %in% Late_stage", RegulonList = upregulated_Late_Germany_symbols$hgnc_symbol)
up.EA.Russia.early <- TCGAanalyze_EAcomplete(TFname = "Russia %in% Early_stage", RegulonList = upregulated_Early_Russia_symbols$hgnc_symbol)
up.EA.Russia.late <- TCGAanalyze_EAcomplete(TFname = "Russia %in% Late_stage", RegulonList = upregulated_Late_Russia_symbols$hgnc_symbol)
up.EA.Poland.early <- TCGAanalyze_EAcomplete(TFname = "Poland %in% Early_stage", RegulonList = upregulated_Early_Poland_symbols$hgnc_symbol)
up.EA.Poland.late <- TCGAanalyze_EAcomplete(TFname = "Poland %in% Late_stage", RegulonList = upregulated_Late_Poland_symbols$hgnc_symbol)
up.EA.Ukraine.early <- TCGAanalyze_EAcomplete(TFname = "Ukraine %in% Early_stage", RegulonList = upregulated_Early_Ukraine_symbols$hgnc_symbol)
up.EA.Ukraine.late <- TCGAanalyze_EAcomplete(TFname = "Ukraine %in% Late_stage", RegulonList = upregulated_Late_Ukraine_symbols$hgnc_symbol)
up.EA.Vietnam.early <- TCGAanalyze_EAcomplete(TFname = "Vietnam %in% Early_stage", RegulonList = upregulated_Early_Vietnam_symbols$hgnc_symbol)
up.EA.Vietnam.late <- TCGAanalyze_EAcomplete(TFname = "Vietnam %in% Late_stage", RegulonList = upregulated_Late_Vietnam_symbols$hgnc_symbol)
up.EA.Korea_South.early <- TCGAanalyze_EAcomplete(TFname = "Korea_South %in% Early_stage", RegulonList = upregulated_Early_Korea_South_symbols$hgnc_symbol)
up.EA.Korea_South.late <- TCGAanalyze_EAcomplete(TFname = "Korea_South %in% Late_stage", RegulonList = upregulated_Late_Korea_South_symbols$hgnc_symbol)

down.EA.US.early <- TCGAanalyze_EAcomplete(TFname = "United_States %in% Early_stage", RegulonList = downregulated_Early_US_symbols$hgnc_symbol)
down.EA.US.late <- TCGAanalyze_EAcomplete(TFname = "United_States %in% Late_stage", RegulonList = downregulated_Late_US_symbols$hgnc_symbol)
down.EA.Germany.early <- TCGAanalyze_EAcomplete(TFname = "Germany %in% Early_stage", RegulonList = downregulated_Early_Germany_symbols$hgnc_symbol)
down.EA.Germany.late <- TCGAanalyze_EAcomplete(TFname = "Germany %in% Late_stage", RegulonList = downregulated_Late_Germany_symbols$hgnc_symbol)
down.EA.Russia.early <- TCGAanalyze_EAcomplete(TFname = "Russia %in% Early_stage", RegulonList = downregulated_Early_Russia_symbols$hgnc_symbol)
down.EA.Russia.late <- TCGAanalyze_EAcomplete(TFname = "Russia %in% Late_stage", RegulonList = downregulated_Late_Russia_symbols$hgnc_symbol)
down.EA.Poland.early <- TCGAanalyze_EAcomplete(TFname = "Poland %in% Early_stage", RegulonList = downregulated_Early_Poland_symbols$hgnc_symbol)
down.EA.Poland.late <- TCGAanalyze_EAcomplete(TFname = "Poland %in% Late_stage", RegulonList = downregulated_Late_Poland_symbols$hgnc_symbol)
down.EA.Ukraine.early <- TCGAanalyze_EAcomplete(TFname = "Ukraine %in% Early_stage", RegulonList = downregulated_Early_Ukraine_symbols$hgnc_symbol)
down.EA.Ukraine.late <- TCGAanalyze_EAcomplete(TFname = "Ukraine %in% Late_stage", RegulonList = downregulated_Late_Ukraine_symbols$hgnc_symbol)
down.EA.Vietnam.early <- TCGAanalyze_EAcomplete(TFname = "Vietnam %in% Early_stage", RegulonList = downregulated_Early_Vietnam_symbols$hgnc_symbol)
down.EA.Vietnam.late <- TCGAanalyze_EAcomplete(TFname = "Vietnam %in% Late_stage", RegulonList = downregulated_Late_Vietnam_symbols$hgnc_symbol)
down.EA.Korea_South.early <- TCGAanalyze_EAcomplete(TFname = "Korea_South %in% Early_stage", RegulonList = downregulated_Early_Korea_South_symbols$hgnc_symbol)
down.EA.Korea_South.late <- TCGAanalyze_EAcomplete(TFname = "Korea_South %in% Late_stage", RegulonList = downregulated_Late_Korea_South_symbols$hgnc_symbol)

#Visualization of the enrichment analysis
# Upregulated Genes 
up.US.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.US.early$ResBP),
                                       GOBPTab = up.EA.US.early$ResBP,
                                       GOCCTab = up.EA.US.early$ResCC,
                                       GOMFTab = up.EA.US.early$ResMF,
                                       PathTab = up.EA.US.early$ResPat,
                                       nRGTab = up.EA.US.early,
                                       nBar = 20,
                                       text.size = 2,
                                       fig.width = 30,
                                       fig.height = 15)

up.US.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.US.late$ResBP),
                                      GOBPTab = up.EA.US.late$ResBP,
                                      GOCCTab = up.EA.US.late$ResCC,
                                      GOMFTab = up.EA.US.late$ResMF,
                                      PathTab = up.EA.US.late$ResPat,
                                      nRGTab = up.EA.US.late,
                                      nBar = 20,
                                      text.size = 2,
                                      fig.width = 30,
                                      fig.height = 15)

up.Germany.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Germany.early$ResBP),
                                            GOBPTab = up.EA.Germany.early$ResBP,
                                            GOCCTab = up.EA.Germany.early$ResCC,
                                            GOMFTab = up.EA.Germany.early$ResMF,
                                            PathTab = up.EA.Germany.early$ResPat,
                                            nRGTab = up.EA.Germany.early,
                                            nBar = 20,
                                            text.size = 2,
                                            fig.width = 30,
                                            fig.height = 15)

up.Germany.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Germany.late$ResBP),
                                           GOBPTab = up.EA.Germany.late$ResBP,
                                           GOCCTab = up.EA.Germany.late$ResCC,
                                           GOMFTab = up.EA.Germany.late$ResMF,
                                           PathTab = up.EA.Germany.late$ResPat,
                                           nRGTab = up.EA.Germany.late,
                                           nBar = 20,
                                           text.size = 2,
                                           fig.width = 30,
                                           fig.height = 15)

up.Russia.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Russia.early$ResBP),
                                           GOBPTab = up.EA.Russia.early$ResBP,
                                           GOCCTab = up.EA.Russia.early$ResCC,
                                           GOMFTab = up.EA.Russia.early$ResMF,
                                           PathTab = up.EA.Russia.early$ResPat,
                                           nRGTab = up.EA.Russia.early,
                                           nBar = 20,
                                           text.size = 2,
                                           fig.width = 30,
                                           fig.height = 15)

up.Russia.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Russia.late$ResBP),
                                          GOBPTab = up.EA.Russia.late$ResBP,
                                          GOCCTab = up.EA.Russia.late$ResCC,
                                          GOMFTab = up.EA.Russia.late$ResMF,
                                          PathTab = up.EA.Russia.late$ResPat,
                                          nRGTab = up.EA.Russia.late,
                                          nBar = 20,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)

up.Poland.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Poland.early$ResBP),
                                           GOBPTab = up.EA.Poland.early$ResBP,
                                           GOCCTab = up.EA.Poland.early$ResCC,
                                           GOMFTab = up.EA.Poland.early$ResMF,
                                           PathTab = up.EA.Poland.early$ResPat,
                                           nRGTab = up.EA.Poland.early,
                                           nBar = 20,
                                           text.size = 2,
                                           fig.width = 30,
                                           fig.height = 15)

up.Poland.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Poland.late$ResBP),
                                          GOBPTab = up.EA.Poland.late$ResBP,
                                          GOCCTab = up.EA.Poland.late$ResCC,
                                          GOMFTab = up.EA.Poland.late$ResMF,
                                          PathTab = up.EA.Poland.late$ResPat,
                                          nRGTab = up.EA.Poland.late,
                                          nBar = 20,
                                          text.size = 2,
                                          fig.width = 30,
                                          fig.height = 15)
up.Ukraine.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Ukraine.early$ResBP),
                                            GOBPTab = up.EA.Ukraine.early$ResBP,
                                            GOCCTab = up.EA.Ukraine.early$ResCC,
                                            GOMFTab = up.EA.Ukraine.early$ResMF,
                                            PathTab = up.EA.Ukraine.early$ResPat,
                                            nRGTab = up.EA.Ukraine.early,
                                            nBar = 20,
                                            text.size = 2,
                                            fig.width = 30,
                                            fig.height = 15)

up.Ukraine.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Ukraine.late$ResBP),
                                           GOBPTab = up.EA.Ukraine.late$ResBP,
                                           GOCCTab = up.EA.Ukraine.late$ResCC,
                                           GOMFTab = up.EA.Ukraine.late$ResMF,
                                           PathTab = up.EA.Ukraine.late$ResPat,
                                           nRGTab = up.EA.Ukraine.late,
                                           nBar = 20,
                                           text.size = 2,
                                           fig.width = 30,
                                           fig.height = 15)

up.Vietnam.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Vietnam.early$ResBP),
                                            GOBPTab = up.EA.Vietnam.early$ResBP,
                                            GOCCTab = up.EA.Vietnam.early$ResCC,
                                            GOMFTab = up.EA.Vietnam.early$ResMF,
                                            PathTab = up.EA.Vietnam.early$ResPat,
                                            nRGTab = up.EA.Vietnam.early,
                                            nBar = 20,
                                            text.size = 2,
                                            fig.width = 30,
                                            fig.height = 15)

up.Vietnam.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Vietnam.late$ResBP),
                                           GOBPTab = up.EA.Vietnam.late$ResBP,
                                           GOCCTab = up.EA.Vietnam.late$ResCC,
                                           GOMFTab = up.EA.Vietnam.late$ResMF,
                                           PathTab = up.EA.Vietnam.late$ResPat,
                                           nRGTab = up.EA.Vietnam.late,
                                           nBar = 20,
                                           text.size = 2,
                                           fig.width = 30,
                                           fig.height = 15)

up.Korea_South.early <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Korea_South.early$ResBP),
                                                GOBPTab = up.EA.Korea_South.early$ResBP,
                                                GOCCTab = up.EA.Korea_South.early$ResCC,
                                                GOMFTab = up.EA.Korea_South.early$ResMF,
                                                PathTab = up.EA.Korea_South.early$ResPat,
                                                nRGTab = up.EA.Korea_South.early,
                                                nBar = 20,
                                                text.size = 2,
                                                fig.width = 30,
                                                fig.height = 15)

up.Korea_South.late <- TCGAvisualize_EAbarplot(tf = rownames(up.EA.Korea_South.late$ResBP),
                                               GOBPTab = up.EA.Korea_South.late$ResBP,
                                               GOCCTab = up.EA.Korea_South.late$ResCC,
                                               GOMFTab = up.EA.Korea_South.late$ResMF,
                                               PathTab = up.EA.Korea_South.late$ResPat,
                                               nRGTab = up.EA.Korea_South.late,
                                               nBar = 20,
                                               text.size = 2,
                                               fig.width = 30,
                                               fig.height = 15)

# Downregulated Genes 
down.US.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.US.early$ResBP),
                                         GOBPTab = down.EA.US.early$ResBP,
                                         GOCCTab = down.EA.US.early$ResCC,
                                         GOMFTab = down.EA.US.early$ResMF,
                                         PathTab = down.EA.US.early$ResPat,
                                         nRGTab = down.EA.US.early,
                                         nBar = 20,
                                         text.size = 2,
                                         fig.width = 30,
                                         fig.height = 15)

down.US.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.US.late$ResBP),
                                        GOBPTab = down.EA.US.late$ResBP,
                                        GOCCTab = down.EA.US.late$ResCC,
                                        GOMFTab = down.EA.US.late$ResMF,
                                        PathTab = down.EA.US.late$ResPat,
                                        nRGTab = down.EA.US.late,
                                        nBar = 20,
                                        text.size = 2,
                                        fig.width = 30,
                                        fig.height = 15)

down.Germany.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Germany.early$ResBP),
                                              GOBPTab = down.EA.Germany.early$ResBP,
                                              GOCCTab = down.EA.Germany.early$ResCC,
                                              GOMFTab = down.EA.Germany.early$ResMF,
                                              PathTab = down.EA.Germany.early$ResPat,
                                              nRGTab = down.EA.Germany.early,
                                              nBar = 20,
                                              text.size = 2,
                                              fig.width = 30,
                                              fig.height = 15)

down.Germany.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Germany.late$ResBP),
                                             GOBPTab = down.EA.Germany.late$ResBP,
                                             GOCCTab = down.EA.Germany.late$ResCC,
                                             GOMFTab = down.EA.Germany.late$ResMF,
                                             PathTab = down.EA.Germany.late$ResPat,
                                             nRGTab = down.EA.Germany.late,
                                             nBar = 20,
                                             text.size = 2,
                                             fig.width = 30,
                                             fig.height = 15)

down.Russia.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Russia.early$ResBP),
                                             GOBPTab = down.EA.Russia.early$ResBP,
                                             GOCCTab = down.EA.Russia.early$ResCC,
                                             GOMFTab = down.EA.Russia.early$ResMF,
                                             PathTab = down.EA.Russia.early$ResPat,
                                             nRGTab = down.EA.Russia.early,
                                             nBar = 20,
                                             text.size = 2,
                                             fig.width = 30,
                                             fig.height = 15)

down.Russia.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Russia.late$ResBP),
                                            GOBPTab = down.EA.Russia.late$ResBP,
                                            GOCCTab = down.EA.Russia.late$ResCC,
                                            GOMFTab = down.EA.Russia.late$ResMF,
                                            PathTab = down.EA.Russia.late$ResPat,
                                            nRGTab = down.EA.Russia.late,
                                            nBar = 20,
                                            text.size = 2,
                                            fig.width = 30,
                                            fig.height = 15)

down.Poland.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Poland.early$ResBP),
                                             GOBPTab = down.EA.Poland.early$ResBP,
                                             GOCCTab = down.EA.Poland.early$ResCC,
                                             GOMFTab = down.EA.Poland.early$ResMF,
                                             PathTab = down.EA.Poland.early$ResPat,
                                             nRGTab = down.EA.Poland.early,
                                             nBar = 20,
                                             text.size = 2,
                                             fig.width = 30,
                                             fig.height = 15)

down.Poland.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Poland.late$ResBP),
                                            GOBPTab = down.EA.Poland.late$ResBP,
                                            GOCCTab = down.EA.Poland.late$ResCC,
                                            GOMFTab = down.EA.Poland.late$ResMF,
                                            PathTab = down.EA.Poland.late$ResPat,
                                            nRGTab = down.EA.Poland.late,
                                            nBar = 20,
                                            text.size = 2,
                                            fig.width = 30,
                                            fig.height = 15)

down.Ukraine.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Ukraine.early$ResBP),
                                              GOBPTab = down.EA.Ukraine.early$ResBP,
                                              GOCCTab = down.EA.Ukraine.early$ResCC,
                                              GOMFTab = down.EA.Ukraine.early$ResMF,
                                              PathTab = down.EA.Ukraine.early$ResPat,
                                              nRGTab = down.EA.Ukraine.early,
                                              nBar = 20,
                                              text.size = 2,
                                              fig.width = 30,
                                              fig.height = 15)

down.Ukraine.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Ukraine.late$ResBP),
                                             GOBPTab = down.EA.Ukraine.late$ResBP,
                                             GOCCTab = down.EA.Ukraine.late$ResCC,
                                             GOMFTab = down.EA.Ukraine.late$ResMF,
                                             PathTab = down.EA.Ukraine.late$ResPat,
                                             nRGTab = down.EA.Ukraine.late,
                                             nBar = 20,
                                             text.size = 2,
                                             fig.width = 30,
                                             fig.height = 15)

down.Vietnam.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Vietnam.early$ResBP),
                                              GOBPTab = down.EA.Vietnam.early$ResBP,
                                              GOCCTab = down.EA.Vietnam.early$ResCC,
                                              GOMFTab = down.EA.Vietnam.early$ResMF,
                                              PathTab = down.EA.Vietnam.early$ResPat,
                                              nRGTab = down.EA.Vietnam.early,
                                              nBar = 20,
                                              text.size = 2,
                                              fig.width = 30,
                                              fig.height = 15)

down.Vietnam.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Vietnam.late$ResBP),
                                             GOBPTab = down.EA.Vietnam.late$ResBP,
                                             GOCCTab = down.EA.Vietnam.late$ResCC,
                                             GOMFTab = down.EA.Vietnam.late$ResMF,
                                             PathTab = down.EA.Vietnam.late$ResPat,
                                             nRGTab = down.EA.Vietnam.late,
                                             nBar = 20,
                                             text.size = 2,
                                             fig.width = 30,
                                             fig.height = 15)

down.Korea_South.early <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Korea_South.early$ResBP),
                                                  GOBPTab = down.EA.Korea_South.early$ResBP,
                                                  GOCCTab = down.EA.Korea_South.early$ResCC,
                                                  GOMFTab = down.EA.Korea_South.early$ResMF,
                                                  PathTab = down.EA.Korea_South.early$ResPat,
                                                  nRGTab = down.EA.Korea_South.early,
                                                  nBar = 20,
                                                  text.size = 2,
                                                  fig.width = 30,
                                                  fig.height = 15)

down.Korea_South.late <- TCGAvisualize_EAbarplot(tf = rownames(down.EA.Korea_South.late$ResBP),
                                                 GOBPTab = down.EA.Korea_South.late$ResBP,
                                                 GOCCTab = down.EA.Korea_South.late$ResCC,
                                                 GOMFTab = down.EA.Korea_South.late$ResMF,
                                                 PathTab = down.EA.Korea_South.late$ResPat,
                                                 nRGTab = down.EA.Korea_South.late,
                                                 nBar = 20,
                                                 text.size = 2,
                                                 fig.width = 30,
                                                 fig.height = 15)

#Save as CSV for GO Analysis
write.csv(upregulated_Early_US_symbols, file = "upregulated_Early_US_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_US_symbols, file = "downregulated_Early_US_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_US_symbols, file = "upregulated_Late_US_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_US_symbols, file = "downregulated_Late_US_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Germany_symbols, file = "upregulated_Early_Germany_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Germany_symbols, file = "downregulated_Early_Germany_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Germany_symbols, file = "upregulated_Late_Germany_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Germany_symbols, file = "downregulated_Late_Germany_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Russia_symbols, file = "upregulated_Early_Russia_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Russia_symbols, file = "downregulated_Early_Russia_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Russia_symbols, file = "upregulated_Late_Russia_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Russia_symbols, file = "downregulated_Late_Russia_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Poland_symbols, file = "upregulated_Early_Poland_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Poland_symbols, file = "downregulated_Early_Poland_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Poland_symbols, file = "upregulated_Late_Poland_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Poland_symbols, file = "downregulated_Late_Poland_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Ukraine_symbols, file = "upregulated_Early_Ukraine_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Ukraine_symbols, file = "downregulated_Early_Ukraine_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Ukraine_symbols, file = "upregulated_Late_Ukraine_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Ukraine_symbols, file = "downregulated_Late_Ukraine_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Vietnam_symbols, file = "upregulated_Early_Vietnam_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Vietnam_symbols, file = "downregulated_Early_Vietnam_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Vietnam_symbols, file = "upregulated_Late_Vietnam_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Vietnam_symbols, file = "downregulated_Late_Vietnam_symbols.csv", row.names = FALSE)

write.csv(upregulated_Early_Korea_South_symbols, file = "upregulated_Early_Korea_South_symbols.csv", row.names = FALSE)
write.csv(downregulated_Early_Korea_South_symbols, file = "downregulated_Early_Korea_South_symbols.csv", row.names = FALSE)

write.csv(upregulated_Late_Korea_South_symbols, file = "upregulated_Late_Korea_South_symbols.csv", row.names = FALSE)
write.csv(downregulated_Late_Korea_South_symbols, file = "downregulated_Late_Korea_South_symbols.csv", row.names = FALSE)

