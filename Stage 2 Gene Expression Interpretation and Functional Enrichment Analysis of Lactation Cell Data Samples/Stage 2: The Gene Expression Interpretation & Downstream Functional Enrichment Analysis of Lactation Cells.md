# **The Gene Expression Interpretation & Downstream Functional Enrichment Analysis of Lactation Cells**

Authors: Nour Nahtay (@NourNahtay), Merylin Ogunlola (@MerylinO) Princess Beatrice Sunday-Jimmy (@LaidyCharm), Onare Opeyemi Mary (@Onare)

The aim of this paper is to interpret a dataset of 590+ gene expression results of lactation cells under different conditions.

Dataset: https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv

## **Heatmap Generation**

We started off by loading and prepping the sample to generate a heatmap. We used two different color palettes, which gave the following results: 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXedKXuYM-hEneadDyKc0Zx6rVI_SPfTEOK1NWAXvPaKV3FghnFmkGyRmYkKWv0iPShka8CrcU5r5O_oEROaaKtcW3sDthiprz4sosdzzCaFVAKQV4QPzCvzddQ2UfFgcP7sRYb-SNBQFcjGmKfnc3J3o0PB?key=DtMktn_nuLAIngoboGJAwA)

_Figure 1: Heatmap with Diverging Palette_

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfIuQfSsEHRsM6lh8SWhvwH7PrMLHqS8aKDY4yhLBQ6K1JsVQkvw5D3J9sxI2VECHycvD0cBKiD9ehaLfJvRdS6eQqou4NK0kGjFjBgnqu2-EuhSH0eScvs1PeryXQ5SYGUrMezStyTzcfHHLGtib7uOigh?key=DtMktn_nuLAIngoboGJAwA)

_Figure 2: Heatmap with Sequential Palette_

The heatmap with the diverging palette offers a more pleasant-to-the-eye interface that allows the viewer to easily distinguish between the over and underexpressed genes. Although the differing hues of the sequential palette could somewhat give insight on the changes, operating with one color is not preferable. 


## **Functional Enrichment Analysis**

To prepare the data for functional enrichment analysis, we began by splitting the data into two groups, then retrieving the mean for each data in order to calculate the fold changes. We then calculated the P-values and used them - along with the fold changes - to subset the upregulated and down regulated genes. 

Upregulated: (fold change > 2 & P-Values < 0.2) & Downregulated: (fold change < -2 & p values < 0.2)

For this analysis, we decided to focus on the upregulated genes, and to perform functional enrichment analysis using PANTHER _(Mi & Thomas)_, then ShinyGo _(Ge, Jung & Yao)_ for confirmation. We used the GO database and an FDR of <0.05 for both websites to ensure consistency and high significance.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdawvKFZC0icpoxg237ekTDEeexyWiAZ3DNaiuqhhWcnsaixFNF2mTsUE6OuNNKJc0CUCq9Oq3jUwUwetL8FLGEu_Fvgjha1abdQ0wmmQ1AN4i4110aRu9oizGqjkWxhSf_ENKHjJlLv2fchFyYRgbh_4b1?key=DtMktn_nuLAIngoboGJAwA)

_Table 1: Functional Enrichment Analysis Using PANTHER_

__![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcoAUfhjEb_VoF_voqfb2GDrbSZXK1dtIaTFQ_xB1BRVBDdSXZphSXAabuNLN9BdEBRevO0JJFK4GNv4DjZftchn2ENzRK1FIt1cCREHeN8RFrB8TX6RB3NQ63LUCwaM8OqwigrzG25WqbRoFW_CTp0HyQ?key=DtMktn_nuLAIngoboGJAwA)__

_Table 2: Functional Enrichment Analysis Using ShinyGO_

The top 5 pathways, number of genes involved, and the negative log 10 of FDR were extracted and used for visualization in the form of a bubble plot

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXckKL78-IOFmOOXPZswBM64uHr6sVcgmFKbWm6u2u_By2gul1V9e0DlCLX_W_0wkEomzEiDymxLLYOY_n_qNtemHNjjD6nw1F-_2F8irI9yhvKuUfSZ9gbx9v8DXSPGzAL4JZze-wNFy_LBycsi_NM79iwd?key=DtMktn_nuLAIngoboGJAwA)

_Table 3: Visualization Data_

__![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfcFV9oRt1kAuoylyG5VwBpeoCEz7aMgI1wlmxcfJ22a8Q8rnH4EXJegDtcOs31s1yPEeu_Qfld3hngg5LDtGbA8ih4SgLcgsNy1Rrh4Thlk1kDLLLbsJuIM_O0jPeQo4rjhBqecIDnXUKlpb0MCxS_7ajZ?key=DtMktn_nuLAIngoboGJAwA)__

_Figure 3: Bubble Plot of the Visualization Data_


## **The Top 3 Pathways**
### Pathway 1: RNA Metabolism Regulation
- **Description:** RNA metabolism regulation encompasses the processes that manage the synthesis, processing, modification, and degradation of RNA molecules, which are vital for protein production and gene expression. This pathway is essential for ensuring RNA levels are properly balanced in cells so that proteins are synthesized at appropriate times and quantities (Wilusz & Sharp, 2013).

- **Relevance:** In this study, genes associated with RNA metabolism were notably 
enriched. This suggests that irregularities in RNA metabolism regulation may play a 
significant role in the condition under investigation, such as cancer or 
developmental disorders.
- **Significance:** With an FDR of 0.0193, this pathway is considered highly 
significant. The enrichment observed might indicate that disruptions in RNA 
metabolism are a key factor driving the observed changes in gene expression, 
potentially highlighting targets for therapeutic strategies



#### **Citations:**

1. Mi H, Thomas P. PANTHER pathway: an ontology-based pathway database coupled with data analysis tools. Methods Mol Biol. 2009;563:123-40. doi: 10.1007/978-1-60761-175-2\_7. PMID: 19597783; PMCID: PMC6608593.

2. Ge SX, Jung D, Yao R. ShinyGO: a graphical gene-set enrichment tool for animals and plants. Bioinformatics. 2020 Apr 15;36(8):2628-2629. doi: 10.1093/bioinformatics/btz931. PMID: 31882993; PMCID: PMC7178415.

#### **Word Count:**
