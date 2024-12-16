# **IDH Status-Based Comparative Gene Expression and Functional Enrichment in Diffuse Glioma Samples**
 
## **Description**
This project explores the role that IDH mutation plays in Glioma. It is split into multiple parts
1. Differential Expression Analysis
2. Functional enrichment analysis
3. k-NN Analysis
4. Random Forest

Dataset:TCGA-LGG

*Biomarker Discovery*

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
2. Treatment level DEA is performed
3. Ouput is used as data to generate heatmaps comparing the clustering of males vs females in both tumor types

### **Functional Enrichment Analysis**
1. The up and downregulated genes for all 4 subsets are extracted individually
2. The data undergoes functional enrichment analysis (particularly gene ontology), selecting the top 5 pathways.
3. Results are visualized in a barplot form and the PDF is saved on the computed

*Machine Learning*
1. Perform k-NN analysis
2. Perform Random Forest

Happy Analyzing!
