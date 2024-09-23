#Install the latest version of BiocManager, which helps install bioconductor tools
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install the latest version
BiocManager::install(version = "3.18") # Install the latest version of Bioconductor
library("BiocManager")

#Install TCGAbiolinks, a package designed to access, analyze, and visualize data from The Cancer Genome Atlas (TCGA) project.
BiocManager::install("TCGAbiolinks")
library("TCGAbiolinks")

#Install limma & edgeR commonly used for analyzing gene expression data, especially for RNA-seq and microarray data, to find differentially expressed genes (DEGs).
BiocManager::install("limma")
BiocManager::install("edgeR")
library("limma")
library("edgeR")

#Install EDASeq for the preprocessing and normalization of RNA-seq data.
BiocManager::install("EDASeq")
library("EDASeq")

#Install gplots for plotting graphs
install.packages('gplots')
library("gplots")

#Install sesameData, which provides preprocessed and annotation data to be used with the sesame package, which is designed for analyzing DNA methylation data, especially from Illumina methylation arrays (e.g., EPIC and 450K).
BiocManager::install("sesameData")
library("sesameData")

# Install SummarizedExperiment, which provides a standardized data structure for storing and managing high-throughput genomic data
BiocManager::install("SummarizedExperiment")
library("SummarizedExperiment")

#Get an Overview of the data on Melanoma
getProjectSummary("TCGA-SKCM")

#Check the other datasets available to work on
?GDCquery

#Query the Data
SKCM <- GDCquery(project = "TCGA-SKCM",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

#Download the dataset
GDCdownload(SKCM)

#Prepare the dataset
SKCM.Data <- GDCprepare(SKCM)

#Check the output of the data
head(SKCM.Data)
View(SKCM.Data)

#Explore the data
SKCM.Data$gender
SKCM.Data$tumor_descriptor
SKCM.Data$ #Insert anything after the dollar sign to check a subset you want, like race, gender, etc
  
#Check the Subsets of the metadata and add them in one data frame
SKCMMETA <- data.frame("gender" = SKCM.Data$gender,
                       "tumor_type" = SKCM.Data$tumor_descriptor,
                       "Barcode" = SKCM.Data$barcode)

#Extract the raw data from the prepared dataset
SKCMRaw <- assays(SKCM.Data)

#Select unstranded data
dim(SKCMRaw$unstranded)
View(SKCMRaw$unstranded)

#Downsize data to 20 primary and metastatic data from both male and female
SelectedBarcodes <- c(subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[c(1:10)])

SelectedData <- SKCMRaw$unstranded[, c(SelectedBarcodes)]
dim(SelectedData)
View(SelectedData) #To check for consistency in expression levels

#Data normalization and filtering
normalized <- TCGAanalyze_Normalization(tabDF = SelectedData, geneInfo = geneInfoHT, method = "geneLength")

#Then filter
filtered <- TCGAanalyze_Filtering(tabDF = normalized,
                                  method = "quantile",
                                  qnt.cut = 0.25)

View(filtered)
dim(filtered)


# Create matrices for pairwise comparison

# For comparison between metastatic data
mat1_MetaMale <- filtered[, subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[1:10]]
mat2_MetaFemale <- filtered[, subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[1:10]]

# Perform DEA for male Metastatic vs female Metastatic
results_metastatic <- TCGAanalyze_DEA(mat1 = mat1_MetaMale,
                                mat2 = mat2_MetaFemale,
                                Cond1type = "male Metastatic",
                                Cond2type = "female Metastatic",
                                pipeline = "edgeR",
                                fdr.cut = 0.01,
                                logFC.cut = 2)

# For comparison between primary
mat1_PrimMale <- filtered[, subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[1:10]]
mat2_PrimFemale <- filtered[, subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[1:10]]

# Perform DEA for male primary vs female Primary
results_primary <- TCGAanalyze_DEA(mat1 = mat1_PrimMale,
                                  mat2 = mat2_PrimFemale,
                                  Cond1type = "male Primary",
                                  Cond2type = "female Primary",
                                  pipeline = "edgeR",
                                  fdr.cut = 0.01,
                                  logFC.cut = 2)

#Volcano Plot
plot(results_metastatic$logFC, -log10(results_metastatic$FDR))
plot(results_primary$logFC, -log10(results_primary$FDR))


#DEA with treatment levels
results_primary.level <- TCGAanalyze_LevelTab(results_primary,"male Primary","female Primary", mat1_PrimMale, mat2_PrimFemale)
results_metastatic.level <- TCGAanalyze_LevelTab(results_metastatic,"male Metastatic","female Metastatic", mat1_MetaMale, mat2_MetaFemale)

head(results_metastatic.level)
head(results_primary.level)

# Create vectors for the barcodes
male_meta_barcodes <- subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[1:10]
female_meta_barcodes <- subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[1:10]
male_prim_barcodes <- subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[1:10]
female_prim_barcodes <- subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[1:10]

# Combine the barcodes for heatmap data
selected_barcodes_meta <- c(male_meta_barcodes, female_meta_barcodes)
selected_barcodes_prim <- c(male_prim_barcodes, female_prim_barcodes)

# Subset the heatmap data
heat.data.meta <- filtered[rownames(results_metastatic.level), selected_barcodes_meta]
heat.data.prim <- filtered[rownames(results_primary.level), selected_barcodes_prim]

# Create a color vector for the samples
sample_colors <- c(rep("royalblue4", 10), rep("thistle2", 10)) # Adjust this if you have different sample counts

#Visualization of DEA
heatmap.2(x = as.matrix(heat.data.meta),
          col = hcl.colors(10, palette = 'Purple-Green'),
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1, #Size of Row and Column
          main = "Heatmap of Metastatic Tumor in Male vs Female",
          na.color = 'black',
          ColSideColors = sample_colors)


#Add a Legend
legend("topright", legend = c("Male", "Female"),
       fill = c("royalblue4", "thistle2"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

heatmap.2(x = as.matrix(heat.data.prim),
          col = hcl.colors(10, palette = 'Purple-Green'),
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1, #Size of Row and Column
          main = "Heatmap of Primary Tumors in Male vs Female",
          na.color = 'black',
          ColSideColors = sample_colors)

#Add a Legend
legend("topright", legend = c("Male", "Female"),
       fill = c("royalblue4", "thistle2"),
       border = "black",  # Black outline
       bty = "n")  # No box around the legend

#Functional Enrichment Analysis

#View volcano plot
plot(results_metastatic$logFC, -log10(results_metastatic$FDR))
plot(results_primary$logFC, -log10(results_primary$FDR))

# Upregulated genes for male and female metastatic
upreg_meta_male <- rownames(subset(results_metastatic.level, logFC > 2 & `male Metastatic` > 0))
upreg_meta_female <- rownames(subset(results_metastatic.level, logFC > 2 & `female Metastatic` > 0))

# Downregulated genes for male and female metastatic
downreg_meta_male <- rownames(subset(results_metastatic.level, logFC < -2 & `male Metastatic` > 0))
downreg_meta_female <- rownames(subset(results_metastatic.level, logFC < -2 & `female Metastatic` > 0))

# Upregulated genes for male and female primary
upreg_prim_male <- rownames(subset(results_primary.level, logFC > 2 & `male Primary` > 0))
upreg_prim_female <- rownames(subset(results_primary.level, logFC > 2 & `female Primary` > 0))

# Downregulated genes for male and female primary
downreg_prim_male <- rownames(subset(results_primary.level, logFC < -2 & `male Primary` > 0))
downreg_prim_female <- rownames(subset(results_primary.level, logFC < -2 & `female Primary` > 0))

#Install and load bioMart, which converts 
BiocManager::install("biomaRt")
library("biomaRt")

#Convert Ensemble ID to gene ID using biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols for male and female metastatic upregulated genes
upreg_meta_male <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = upreg_meta_male,
                         mart = mart)$hgnc_symbol

upreg_meta_female <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = upreg_meta_female,
                           mart = mart)$hgnc_symbol

# Convert Ensemble IDs to gene symbols for male and female metastatic downregulated genes
downreg_meta_male <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = downreg_meta_male,
                           mart = mart)$hgnc_symbol

downreg_meta_female <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                             filters = 'ensembl_gene_id',
                             values = downreg_meta_female,
                             mart = mart)$hgnc_symbol

# Convert Ensemble IDs to gene symbols for male and female primary upregulated genes
upreg_prim_male <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = upreg_prim_male,
                         mart = mart)$hgnc_symbol

upreg_prim_female <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = upreg_prim_female,
                           mart = mart)$hgnc_symbol

# Convert Ensemble IDs to gene symbols for male and female primary downregulated genes
downreg_prim_male <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = downreg_prim_male,
                           mart = mart)$hgnc_symbol

downreg_prim_female <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                             filters = 'ensembl_gene_id',
                             values = downreg_prim_female,
                             mart = mart)$hgnc_symbol

# Functional Enrichment Analysis for Metastatic (Upregulated and Downregulated)
up.EA.meta_male <- TCGAanalyze_EAcomplete(TFname = "Upregulated Metastatic Male", RegulonList = upreg_meta_male)
down.EA.meta_male <- TCGAanalyze_EAcomplete(TFname = "Downregulated Metastatic Male", RegulonList = downreg_meta_male)

up.EA.meta_female <- TCGAanalyze_EAcomplete(TFname = "Upregulated Metastatic Female", RegulonList = upreg_meta_female)
down.EA.meta_female <- TCGAanalyze_EAcomplete(TFname = "Downregulated Metastatic Female", RegulonList = downreg_meta_female)

# Functional Enrichment Analysis for Primary (Upregulated and Downregulated)
up.EA.prim_male <- TCGAanalyze_EAcomplete(TFname = "Upregulated Primary Male", RegulonList = upreg_prim_male)
down.EA.prim_male <- TCGAanalyze_EAcomplete(TFname = "Downregulated Primary Male", RegulonList = downreg_prim_male)

up.EA.prim_female <- TCGAanalyze_EAcomplete(TFname = "Upregulated Primary Female", RegulonList = upreg_prim_female)
down.EA.prim_female <- TCGAanalyze_EAcomplete(TFname = "Downregulated Primary Female", RegulonList = downreg_prim_female)

#Visualization of the enrichment analysis of upregulated Metastatic Female Samples
TCGAvisualize_EAbarplot(tf = rownames(up.EA.meta_female$ResBP),
                        GOBPTab = up.EA.meta_female$ResBP,
                        GOCCTab = up.EA.meta_female$ResCC,
                        GOMFTab = up.EA.meta_female$ResMF,
                        PathTab = up.EA.meta_female$ResPat,
                        nRGTab = upreg_meta_female,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of upregulated Metastatic male Samples
TCGAvisualize_EAbarplot(tf = rownames(up.EA.meta_male$ResBP),
                        GOBPTab = up.EA.meta_male$ResBP,
                        GOCCTab = up.EA.meta_male$ResCC,
                        GOMFTab = up.EA.meta_male$ResMF,
                        PathTab = up.EA.meta_male$ResPat,
                        nRGTab = upreg_meta_male,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of downregulated Metastatic male Samples
TCGAvisualize_EAbarplot(tf = rownames(down.EA.meta_male$ResBP),
                        GOBPTab = down.EA.meta_male$ResBP,
                        GOCCTab = down.EA.meta_male$ResCC,
                        GOMFTab = down.EA.meta_male$ResMF,
                        PathTab = down.EA.meta_male$ResPat,
                        nRGTab = downreg_meta_male,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of downregulated Metastatic female Samples
TCGAvisualize_EAbarplot(tf = rownames(down.EA.meta_female$ResBP),
                        GOBPTab = down.EA.meta_female$ResBP,
                        GOCCTab = down.EA.meta_female$ResCC,
                        GOMFTab = down.EA.meta_female$ResMF,
                        PathTab = down.EA.meta_female$ResPat,
                        nRGTab = downreg_meta_female,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of upregulated Primary Female Samples
TCGAvisualize_EAbarplot(tf = rownames(up.EA.prim_female$ResBP),
                        GOBPTab = up.EA.prim_female$ResBP,
                        GOCCTab = up.EA.prim_female$ResCC,
                        GOMFTab = up.EA.prim_female$ResMF,
                        PathTab = up.EA.prim_female$ResPat,
                        nRGTab = upreg_prim_female,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of upregulated Primary male Samples
TCGAvisualize_EAbarplot(tf = rownames(up.EA.prim_male$ResBP),
                        GOBPTab = up.EA.prim_male$ResBP,
                        GOCCTab = up.EA.prim_male$ResCC,
                        GOMFTab = up.EA.prim_male$ResMF,
                        PathTab = up.EA.prim_male$ResPat,
                        nRGTab = upreg_prim_male,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of downregulated Primary male Samples
TCGAvisualize_EAbarplot(tf = rownames(down.EA.prim_male$ResBP),
                        GOBPTab = down.EA.prim_male$ResBP,
                        GOCCTab = down.EA.prim_male$ResCC,
                        GOMFTab = down.EA.prim_male$ResMF,
                        PathTab = down.EA.prim_male$ResPat,
                        nRGTab = downreg_prim_male,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

#Visualization of the enrichment analysis of downregulated Primary female Samples
TCGAvisualize_EAbarplot(tf = rownames(down.EA.prim_female$ResBP),
                        GOBPTab = down.EA.prim_female$ResBP,
                        GOCCTab = down.EA.prim_female$ResCC,
                        GOMFTab = down.EA.prim_female$ResMF,
                        PathTab = down.EA.prim_female$ResPat,
                        nRGTab = downreg_prim_female,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)
