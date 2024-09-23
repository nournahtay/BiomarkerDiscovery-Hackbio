# **Gender-Based Comparative Gene Expression and Functional Enrichment in Melanoma Samples** 

Authors: Nour Nahtay (@NourNahtay), Merylin Ogunlola (@MerylinO), Igwebuike Oluchukwu Vivian (@igwebuikevee0000), Titilola Shittu (@lola), Princess Beatrice Sunday-Jimmy (@LaidyCharm), Onare Opeyemi Mary (@Onare)

The aim of this analysis is to interpret a glioblastoma gene expression dataset and study their pathway involvement.

Dataset: https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv

## **Heatmap Generation**

Dataset is loaded and prepped for heatmap generation. We used two different color palettes, which gave the following results: 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXedKXuYM-hEneadDyKc0Zx6rVI_SPfTEOK1NWAXvPaKV3FghnFmkGyRmYkKWv0iPShka8CrcU5r5O_oEROaaKtcW3sDthiprz4sosdzzCaFVAKQV4QPzCvzddQ2UfFgcP7sRYb-SNBQFcjGmKfnc3J3o0PB?key=DtMktn_nuLAIngoboGJAwA)

_Figure 1: Heatmap with Diverging Palette_

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfIuQfSsEHRsM6lh8SWhvwH7PrMLHqS8aKDY4yhLBQ6K1JsVQkvw5D3J9sxI2VECHycvD0cBKiD9ehaLfJvRdS6eQqou4NK0kGjFjBgnqu2-EuhSH0eScvs1PeryXQ5SYGUrMezStyTzcfHHLGtib7uOigh?key=DtMktn_nuLAIngoboGJAwA)

_Figure 2: Heatmap with Sequential Palette_

The heatmap with the diverging palette offers a more pleasant-to-the-eye interface that allows the viewer to easily distinguish between the over and underexpressed genes, as opposed to the sequential palette, which was not as insightful.


## **Functional Enrichment Analysis**

To perform functional enrichment analysis, the data is split into two groups to calculate the p-values, and the means needed to calculate the fold change. The p-values and fold changes will be used as cutoff to subset up and down regulated genes. 

Functional enrichment for this task was performed on the upregulated genes only using PANTHER _(Mi & Thomas)_, then ShinyGo _(Ge, Jung & Yao)_ for confirmation. We used the GO database and an FDR of <0.05 for both websites to ensure consistency and high significance.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdawvKFZC0icpoxg237ekTDEeexyWiAZ3DNaiuqhhWcnsaixFNF2mTsUE6OuNNKJc0CUCq9Oq3jUwUwetL8FLGEu_Fvgjha1abdQ0wmmQ1AN4i4110aRu9oizGqjkWxhSf_ENKHjJlLv2fchFyYRgbh_4b1?key=DtMktn_nuLAIngoboGJAwA)

_Table 1: Functional Enrichment Analysis Using PANTHER_

__![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcoAUfhjEb_VoF_voqfb2GDrbSZXK1dtIaTFQ_xB1BRVBDdSXZphSXAabuNLN9BdEBRevO0JJFK4GNv4DjZftchn2ENzRK1FIt1cCREHeN8RFrB8TX6RB3NQ63LUCwaM8OqwigrzG25WqbRoFW_CTp0HyQ?key=DtMktn_nuLAIngoboGJAwA)__

_Table 2: Functional Enrichment Analysis Using ShinyGO_

The top 5 pathways, number of genes involved, and the negative log 10 of FDR were extracted and used for visualization in the form of a bubble plot.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXckKL78-IOFmOOXPZswBM64uHr6sVcgmFKbWm6u2u_By2gul1V9e0DlCLX_W_0wkEomzEiDymxLLYOY_n_qNtemHNjjD6nw1F-_2F8irI9yhvKuUfSZ9gbx9v8DXSPGzAL4JZze-wNFy_LBycsi_NM79iwd?key=DtMktn_nuLAIngoboGJAwA)

_Table 3: Visualization Data_

__![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfcFV9oRt1kAuoylyG5VwBpeoCEz7aMgI1wlmxcfJ22a8Q8rnH4EXJegDtcOs31s1yPEeu_Qfld3hngg5LDtGbA8ih4SgLcgsNy1Rrh4Thlk1kDLLLbsJuIM_O0jPeQo4rjhBqecIDnXUKlpb0MCxS_7ajZ?key=DtMktn_nuLAIngoboGJAwA)__

_Figure 3: Bubble Plot of the Visualization Data_


## **The Top 3 Pathways**

### Pathway 1: RNA Metabolism Regulation
In glioblastoma, RNA metabolism is notably disrupted, impacting tumor growth and treatment resistance. Abnormal RNA processing, including faulty splicing and modifications, results in defective proteins that drive tumor progression. GBM cells often exhibit abnormal alternative splicing, which affects crucial cellular processes and promotes malignancy. Additionally, altered RNA modifications can influence mRNA stability and translation, further contributing to tumor aggressiveness. Non-coding RNAs, such as microRNAs and long non-coding RNAs, are also deregulated, affecting gene expression and chromatin structure. _(Dong & Cui, 2019)_

### Pathway 2: RNA Biosynthesis Regulation
RNA biosynthesis involves transcription of pre-mRNA, which undergoes multiple post-transcriptional modifications, including splicing. Splicing removes introns and joins exons, resulting in mature mRNA that is translated into proteins. Kim et al. highlight aberrant RNA splicing in glioblastoma, focusing on the role of the SON protein, which influences the splicing machinery—an essential component of RNA biosynthesis. _(Kim et al., 2021)_

### Pathway 3: Regulation of Transcription by RNA Polymerase II
RNA Polymerase II transcribes DNA into mRNA in eukaryotic cells, and its regulation ensures genes are expressed appropriately. In glioblastoma, G9a enhances the transcription of the TCF12 gene by activating p-STAT3, which recruits RNA Polymerase II to its target promoters. This boosts the transcription of genes involved in tumor progression and radio-resistance. By modulating chromatin accessibility via histone modification, G9a indirectly controls Pol II's efficiency in transcribing oncogenic genes, contributing to glioblastoma's aggressive nature and linking epigenetic regulation to Pol II-driven gene expression. _(Li et al., 2023)_


#### **Citations:**

1. Mi H, Thomas P. PANTHER pathway: an ontology-based pathway database coupled with data analysis tools. Methods Mol Biol. 2009;563:123-40. doi: 10.1007/978-1-60761-175-2\_7. PMID: 19597783; PMCID: PMC6608593.

2. Ge SX, Jung D, Yao R. ShinyGO: a graphical gene-set enrichment tool for animals and plants. Bioinformatics. 2020 Apr 15;36(8):2628-2629. doi: 10.1093/bioinformatics/btz931. PMID: 31882993; PMCID: PMC7178415.
   
3. Dong, Z., & Cui, H. (2019). Epigenetic modulation of metabolism in glioblastoma. Seminars in Cancer Biology, 57, 45–51.
https://doi.org/10.1016/j.semcancer.2018.09.002

4. Kim, J.-H., Jeong, K., Li, J., Murphy, J. M., Vukadin, L., Stone, J. K., Richard, A., Tran, J., Gillespie, G. Y., Flemington, E. K., Sobol, R. W., Lim, S.-T. S., & Ahn, E.-Y. E. 
(2021). SON drives oncogenic RNA splicing in glioblastoma by regulating PTBP1/PTBP2 switching and RBFOX2 activity. Nature Communications, 12(1), 5551.
https://doi.org/10.1038/s41467-021-25892-x

5. Li, X.-L., Xie, Y., Chen, Y.-L., Zhang, Z.-M., Tao, Y.-F
