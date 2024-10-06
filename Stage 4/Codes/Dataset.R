# Load Important Packages
library("BiocManager")
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("EDASeq")
library("gplots")
library("sesameData")
library("SummarizedExperiment")

# Get an overview on  LGG data
getProjectSummary("TCGA-LGG")

# Query the LGG Data
LGG <- GDCquery(project = "TCGA-LGG",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification")

# Download the dataset
GDCdownload(LGG)

# Prepare the dataset
LGG.Data <- GDCprepare(LGG)

# Check the output of the data for LGG
head(LGG.Data)
View(LGG.Data)

# Explore the data
LGG.Data$gender
LGG.Data$barcode
LGG.Data$patient
LGG$
  
  # Check the Subsets of the metadata and add them in one data frame
  LGG.metaData <- data.frame("IDH_status" = LGG.Data$paper_IDH.status,
                             "Expression_cluster" = LGG.Data$paper_Pan.Glioma.RNA.Expression.Cluster,
                             "Barcode" = LGG.Data$barcode)

# Extract the raw data from the prepared dataset
LGG.Raw <- assays(LGG.Data)

#Select unstranded data
dim(LGG.Raw$unstranded)
View(LGG.Raw$unstranded)

# Subset Meta Data
SelectedBarcodes <- c(subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr1")$Barcode,
                      subset(LGG.metaData, IDH_status == "Mutant" & Expression_cluster == "LGr4")$Barcode,
                      subset(LGG.metaData, IDH_status == "WT" & Expression_cluster == "LGr4")$Barcode)

#Retrieve the unstranded data for the selected barcodes
Selected_LGG_Data <- LGG.Raw$unstranded[, c(SelectedBarcodes)]
dim(Selected_LGG_Data)
View(Selected_LGG_Data)

#Normalize Data
LGG_normalized <- TCGAanalyze_Normalization(tabDF = Selected_LGG_Data, geneInfo = geneInfoHT, method = "geneLength")

# Then Filter
LGG_filtered <- TCGAanalyze_Filtering(tabDF = LGG_normalized,
                                      method = "quantile",
                                      qnt.cut = 0.25)

View(LGG_filtered)
dim(LGG_filtered)

