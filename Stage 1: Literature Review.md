### **Large-scale RNA-Seq Transcriptome Analysis of 4043 Cancers and 548 Normal Tissue Controls across 12 TCGA Cancer Types**

Authors: Nour Nahtay (@nournahtay) Ibrahim (@IbrahimFangary) Yara Haitham (@YaraHaitham) Salma Ismail (@Salmaismail28)

@Salma

Cancer is a complex disease driven by gene expression pattern changes due to the accumulation of mutations or epigenetic modifications. In recent years, RNA-Seq has been proven crucial for analyzing transcriptomes to help uncover molecular mechanisms underlying cancer. In this study, researchers conducted a comprehensive analysis of gene expression alteration in 4043 cancer samples and 548 normal tissue controls sourced from TCGA, covering 12 major cancer types.

 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXe4TbygBsP90AzNLMJBK3pX_1N_NLqUfBfFrR2802s7H9L0rqOQsJp4jACQVoaXVgKSk2LR6J0hp9y2bgkgSXWsOwvfShWS2PfOL6VS1kcIysfRIJFGUL0Ml-nAPtWiCYFcdF6bUdDWHbufOhR7LCE_Kyrn?key=r-gOvmhs1nKRqmGTtCDj6A)

_Table 1: Number of Cancer and Normal Samples of the 12 Cancer Types Used in this Analysis_


### @Nour & @Ibrahim

### **Differential Gene Expression Analysis**

This analysis starts by dividing samples into “primary tumor” & “solid tissue normal” groups. Raw counts were compiled and analyzed using the edgeR bioconductor package to find the differentially expressed genes between the two phenotypic groups.


### **Gene Set Association Analysis**

GSAA starts with RNA-seq data normalization using the DEseq tool, followed by clustering with APCluster. The gene sets were validated using the Pearson correlation coefficient. For gene set association analysis, researchers developed their own software, called GSAASeqSP, which was used to perform gene-level differential analysis (Signal2Noise) and gene set association analysis (Weighted\_KS).


### **Pathway enrichment analysis and disease association analysis**

To understand the biological involvement of the selected gene sets, Researchers carried out three types of pathway enrichment analyses and a disease association analysis using WebGestalt. The three types of pathway enrichment analyses were GO, KEGG, and Pathway Commons analyzes.

![Figure 6](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdZzGfWwyVpb835etErlmD0TnDjyjw3Iukh8HTmJxw8VWN9fzTPZrBmCWEMp69Vle4gxuY74wkAs48cxfFG2_mdPGE5Z7-wi5XXBCOB8dXjmJTlkglLgz4eb09zpAiS_QJrW_TdcdzaL8fGlLy8s5XOqVg?key=r-gOvmhs1nKRqmGTtCDj6A)

_Figure 1: Overview of the pipeline_


### @Nour @Ibrahim & @Yara

### **Results** 

The gene set association analysis identified 46 significant gene sets for 7 cancer types, with seven being cross-cancer gene signatures with significantly altered expression levels in at least four cancer types. 

|             |                                                                  |
| ----------- | ---------------------------------------------------------------- |
| Cluster     | Potential Prognostic Marker                                      |
| CLUSTER241  | MKI67, TOP2A, ERBB2/HER2 and Kinesins                            |
| CLUSTER514  | EZH2, AURKA, BIRC5, TK1, PLK1, RAD51, HMMR, CCNB1, PBK and CDKN3 |
| CLUSTER1011 | BRIP1, KIAA1524/CIP2A and CENPH                                  |
| CLUSTER932  | CCNE1, CCNE2, CDC6, E2F7 and UHRF1                               |
| CLUSTER574  | FOXM1, TYMS, RRM2 and SKA1                                       |
| CLUSTER3137 | SKP2, RRM1 and DNMT1                                             |
| CLUSTER184  | AURKB and UBE2C                                                  |

_Table 1: Key Potential Prognostic Marker in Each Cluster_ 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfb2Fa5VIMBuC4adUKmI5kr1BGrmFl6DCmW1rHtm8t-L2VQ5aSl6qlOzP9Hf7E33IZs_vB23FA5aW2S9XeG0TUxl3geqM1ulvmorWwKeXbJi7YwLszeaUUg4X_rRQuSEUNAAIdaUMZ3UMcq08E1JrA9gm8?key=r-gOvmhs1nKRqmGTtCDj6A)

_Figure 2: The significant gene sets with FDR <0.15 across all cancer types_

To validate these sets, the top two most DE genes from these sets were used to create a 14-gene signature and Leave-One-Out Cross-Validation was used to test the feasibility of these signatures in distinguishing between normal and cancerous tissue samples for these cancers. The predictive accuracy of the LOOCV ranged from 88-99%, which confirmed their relevance.

@ Nour 

**Conclusion**

These findings offer major insight into the role of transcripts and their mutations in instigating tumorigenesis and metastasis. By bridging the gap between gene expression alteration and cancer mechanisms, these signatures can serve as foundational elements in developing diagnostics tools and even treatments strategies. 

**Word Count**: 325

**Article:** https\://www\.nature.com/articles/srep13413

**Video:**

Editor: Ibrahim

Voice Over: Salma & Nour

Script: Salma Nour & Yara
