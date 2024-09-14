## **The Gene Expression Interpretation & Downstream Functional Enrichment Analysis of Lactation Cells**

### **Description**
This project is centered around the analysis of a pregnancy and lactation cells gene expression dataset. The project is split into two parts:
1. Interpretation of the gene expression data set using RStudio
2. Functional enrichment analysis on the upregulated genes
3. Visualization of the functional enrichment analysis results using RStudio

Dataset: https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/pregnancyLactationCells.csv

You will need the following packages on RStudio:
1. gplots - Heatmap generation: 
2. dplyr - Data manipulation
3. ggplot2 - Functional enrichment analysis visualization

These packages can be installed using install.packages()

The following functional enrichment analysis tools can be used for this analysis
1. ShinyGO: http://bioinformatics.sdstate.edu/go/
2. GOrilla: https://cbl-gorilla.cs.technion.ac.il/
3. PANTHER: https://geneontology.org/

#### **Gene Expression Interpretation**
1. The URL is to be pasted directly onto RStudio, then read using read.csv
     You can customize the rows and columns to be included using row.names or col.names
2. Heatmaps are generated using heatmpa.2() from gplots. Please note that the dataset is directly converted to the matrix form using as.matrix() 
     The heatmap is completely customizable, from removing dendrograms, to alternating the scale, to the color palette used and more.
     For other color palette options, use https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
     If you recieve an error related to the margains, simply move make the plot display larger with your pointer.


