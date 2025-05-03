# **Uncovering Country-Specific Gene Expression Patterns in Stomach Cancer**
 
## **Description**
This project highlights the differences in differentially expressed genes of samples from different countries. The pipeline is split into 4 essential steps:
1. Subsetting data according to countries
2. Preproccessing data by normalization and filtration
3. Differential expression analysis
4. Gene function analysis

Dataset:TCGA-STAD

You will need the following packages in RStudio:
![68747470733a2f2f6c68372d72742e676f6f676c6575736572636f6e74656e742e636f6d2f646f63737a2f41445f346e5866376371615537676b4363524e70654232643943662d583859664b6d646164324a5f78346733566d50325557666c75637830554b4434506774497](https://github.com/user-attachments/assets/ce2b1181-4157-4b3f-bec5-3c1f3dff3600)

These packages can be installed using install.packages()

### **Preprocessing**
Dataset is uploaded directly onto RStudio using TCGABiolink, then preproccessed by:
  1. Subsetting into the designated groups
![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcgs3b-jBDin7N-x0ncfJdEDN7ag6Qj85sFUTVom1uV3-xd0_x3JSj0U6ZqhFEsFEJVZsfZqntzMU1sy2XCc5MsZPqkGuHcknlN6YaUDpTO4FBcAH5HYzymZn1LBt3-KLHP47HLfYrojNJUfknIRTaBS50?key=AgyQUxUvRaYsce19bm41Iw)
  2. Normalization
  3. Filteration with a quantile cut off 0.25 

### **Differential Expressionn Analysis**
1. Mutant and WT Samples are grouped together
2. Treatment level DEA is performed with a fdr cut = 0.1 and logFC.cut = 2
3. Ouput is used as data to generate heatmaps comparing the clustering of males vs females in both tumor types

### **Functional Enrichment Analysis**
1. The data undergoes functional enrichment analysis (particularly gene ontology), selecting the top 5 pathways.
2. Results are visualized in a barplot form and the PDF is saved on the computer
