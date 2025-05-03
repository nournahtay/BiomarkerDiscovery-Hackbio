#Get an Overview of the data on Melanoma
getProjectSummary("TCGA-SKCM")

#Check the other datasets available to work on
?GDCquery

#Query the data and extract the needed category
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

#Check the Subsets of the metadata and add them in one data frame (You have to drag your mouse over this one)
SKCMMETA <- data.frame("gender" = SKCM.Data$gender,
                       "tumor_type" = SKCM.Data$tumor_descriptor,
                       "Barcode" = SKCM.Data$barcode)

#Extract the raw data from the prepared dataset
SKCMRaw <- assays(SKCM.Data)

#View unstranded data
dim(SKCMRaw$unstranded)
View(SKCMRaw$unstranded)

#Downsize data to 20 primary and metastatic data from both male and female
SelectedBarcodes <- c(subset(SKCMMETA, tumor_type == "Metastatic" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "male")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Metastatic" & gender == "female")$Barcode[c(1:10)],
                      subset(SKCMMETA, tumor_type == "Primary" & gender == "female")$Barcode[c(1:10)])
SelectedData <- SKCMRaw$unstranded[, c(SelectedBarcodes)]

#View data
dim(SelectedData)
View(SelectedData)

#Normalize data
normalized <- TCGAanalyze_Normalization(tabDF = SelectedData, geneInfo = geneInfoHT, method = "geneLength")

#filter
filtered <- TCGAanalyze_Filtering(tabDF = normalized,
                                  method = "quantile",
                                  qnt.cut = 0.25)
#View Data
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

# For comparison between primary data
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

#View data
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

#Upregulated genes for male and female metastatic
upreg_meta <- rownames(subset(results_metastatic.level, logFC > 2))
downreg_meta <- rownames(subset(results_metastatic.level, logFC < -2))
upreg_prim_male <- rownames(subset(results_primary.level, logFC > 2))
downreg_prim_male <- rownames(subset(results_primary.level, logFC < -2))

#Convert Ensemble ID to gene ID using biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensemble IDs to gene symbols for male and female metastatic upregulated genes
upreg_meta <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = upreg_meta,
                         mart = mart)$hgnc_symbol

downreg_meta <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = downreg_meta,
                           mart = mart)$hgnc_symbol

upreg_prim <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = upreg_prim,
                         mart = mart)$hgnc_symbol

downreg_prim <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                           filters = 'ensembl_gene_id',
                           values = downreg_prim,
                           mart = mart)$hgnc_symbol

# Functional Enrichment Analysis for Metastatic (Upregulated and Downregulated)
up.EA.meta <- TCGAanalyze_EAcomplete(TFname = "Upregulated Metastatic", RegulonList = upreg_meta)
down.EA.meta <- TCGAanalyze_EAcomplete(TFname = "Downregulated Metastatic ", RegulonList = downreg_meta)
up.EA.prim <- TCGAanalyze_EAcomplete(TFname = "Upregulated Primary", RegulonList = upreg_prim)
down.EA.prim <- TCGAanalyze_EAcomplete(TFname = "Downregulated Primary", RegulonList = downreg_prim)


#Visualization of the enrichment analysis of upregulated Metastatic Female Samples
TCGAvisualize_EAbarplot(tf = rownames(up.EA.meta$ResBP),
                        GOBPTab = up.EA.meta$ResBP,
                        GOCCTab = up.EA.meta$ResCC,
                        GOMFTab = up.EA.meta$ResMF,
                        PathTab = up.EA.meta$ResPat,
                        nRGTab = upreg_meta,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

TCGAvisualize_EAbarplot(tf = rownames(down.EA.meta$ResBP),
                        GOBPTab = down.EA.meta$ResBP,
                        GOCCTab = down.EA.meta$ResCC,
                        GOMFTab = down.EA.meta$ResMF,
                        PathTab = down.EA.meta$ResPat,
                        nRGTab = downreg_meta,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

TCGAvisualize_EAbarplot(tf = rownames(up.EA.prim$ResBP),
                        GOBPTab = up.EA.prim$ResBP,
                        GOCCTab = up.EA.prim$ResCC,
                        GOMFTab = up.EA.prim$ResMF,
                        PathTab = up.EA.prim$ResPat,
                        nRGTab = upreg_prim,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

TCGAvisualize_EAbarplot(tf = rownames(down.EA.prim$ResBP),
                        GOBPTab = down.EA.prim$ResBP,
                        GOCCTab = down.EA.prim$ResCC,
                        GOMFTab = down.EA.prim$ResMF,
                        PathTab = down.EA.prim$ResPat,
                        nRGTab = downreg_prim,
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)
