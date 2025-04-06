# **IDH Status-Based Comparative Gene Expression and Functional Enrichment in Diffuse Glioma Samples** 

Authors: Nour Nahtay (@NourNahtay), Merylin Ogunlola (@MerylinO), Igwebuike Oluchukwu Vivian (@igwebuikevee0000), Rutuja Pangare (@rutuja0502), ), Titilola Shittu (@lola)**

**late submission

Adult diffuse gliomas are the most frequent malignant tumors affecting the central nervous system. Survival rates differ significantly based on subtype. Low-grade gliomas have a high 5-year survival rate, whereas high-grade gliomas are associated with much lower 5-year survival rates. _(Molinaro AM et al. 2019)_

Although histopathologic classification is well-established it is prone to significant intra- and inter-observer variability, especially in grade II-III tumors. _(Whitfield BT et al.,2022)_ 

One key identifier is The isocitrate dehydrogenase family, which comprises three isoforms located in the cytoplasm and peroxisomes (IDH1), and mitochondria (IDH2 and IDH3). These enzymes are involved in cellular processes including mitochondrial oxidative phosphorylation, glutamine metabolism, and more. _(Cairns R et al.,2013)_

Research suggested classifying gliomas into three groups: 

1. IDH wild-type

2. IDH mutant-codel

3. IDH mutant-non-codel

We sought to determine pan-glioma expression subtypes through unsupervised clustering analysis of  RNA-seq profiles by clustering Gliomas Identifies into four clusters labeled LGr1–4. _(Ceccarelli et al.,2016)_

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfjciKp4uGbSg46hcDPrlCI15K2dMQQe-NqAhK_anJ43sUvjdXj1GV7iaR1IMSZtgXKlp6hM2GKsEsr8o4V6B5v4Dy-Yixf7tjnGhM7skrf15-aiSTOos0bqf3m8dhszyrwehGhmAtIxITIM63o2e5GjJ3K?key=AgyQUxUvRaYsce19bm41Iw)****

Figure 1: An Overview of the Packages Used for Biomarker Discovery Pipeline

Created in BioRender.com


**Data Preprocessing and Analysis**

TCGA-LGG data is uploaded onto RStudio, then filtered according to IDH status and type. For this analysis we decided to compare between LGr1 & LGr2, and LGr1 & LGr4. Normalization and filtering were performed with a cutoff of 0.25.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcgs3b-jBDin7N-x0ncfJdEDN7ag6Qj85sFUTVom1uV3-xd0_x3JSj0U6ZqhFEsFEJVZsfZqntzMU1sy2XCc5MsZPqkGuHcknlN6YaUDpTO4FBcAH5HYzymZn1LBt3-KLHP47HLfYrojNJUfknIRTaBS50?key=AgyQUxUvRaYsce19bm41Iw)

Figure 2: Data subset for this Analysis

**Differential Expression Analysis**

We performed DEA using edgeR to identify DEGs in glioma samples with significance criteria set at an adjusted log2 fold change (log FC) > 2 and a false discovery rate (FDR) < 0.01.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcC7Zw8-etKj6-AFDF71s5c5Fc3ySC7d97AHpOe1OUzy8TmrJI4b8I7g-tziH4fxsFuyaNmkJd5iHOcveETGpes7gSwp_esL0wVraDNHAjSP48ym6SeLBsbUCw67ss59itFDGJqR0bkpmQvoLgc2Vo_BVlc?key=AgyQUxUvRaYsce19bm41Iw)

Figure 3: Heatmap showcasing the DEA results of mutant type glioma in LGr 1 and 4

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcGnBptjSb91g4XWpkQxtsXZzFpF94BOFBqXcRhQkCxx8eeTJlpHFtHZF31RWVwrkxQw28XiISBRXt06JyelXV2vz7JvYyuVnTh-8k4f98FiZ3SmyxNBXuakypz0yFD_Uc5HM5vS17Kd6WKKUcy7p2hcezo?key=AgyQUxUvRaYsce19bm41Iw)

Figure 4: Heatmap showcasing the DEA results of WT glioma in LGr 1 and 4

**Functional Enrichment Analysis**

We extracted the DEGs for each subset then used biomaRt to convert EnsembleIDs to GeneIDS to perform FEA, and more specifically GO analysis to pinpoint the biological processes, molecular functions, and cellular components linked to these genes.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcLoc6o6zRxiOG-Y3HzLi1TSNQFCpDJD9l-9-EssVE_UAcOFaE3SyJ42M-WvEeLpTBLLujNpK5ud9EdRnYpa5PschqaM27CCWAas1dZ62wdCPkie-Gu1ubN6Do8HgLeJnbIRHNIh6zZYLHSTn6gHmrVaiRB?key=AgyQUxUvRaYsce19bm41Iw)

Figure 5: FEA results of the upregulated genes of WT glioma on LGr1

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcPrmajoN-2l1DjKOfEnHehkwU5IEjbENH0X6xF5SMTpIr_ZQQhBSg0td8ywg1mBlucYjG8G0vTCuu-i5bP_U0xBDopFtpHsTQ02OJu5Y2SshAJcgGM0kzSpkloN0WULqUDx_nUOm017iXRmKJ8aXst6g02?key=AgyQUxUvRaYsce19bm41Iw)

Figure 6: FEA results of the upregulated genes of WT glioma on LGr4

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXemaXd6tvyq-Xl836Il1sSvHDChFKGzT3UmVQWxNLwIfb9VM0CmAh2g6qbHqs_dyi6TloZ0GJpi6xo8-c_W7p-cHjQRglj6_Uj9XtSN8LqvI9vOq_PXLV5dL0i2NHPh_qab63bqonI66ggD67SlWvIrRWY?key=AgyQUxUvRaYsce19bm41Iw)

Figure 7: FEA results of the upregulated genes of mutant glioma on LGr1

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeJFLcpvorGtNpOaTMMmgt13c4WsORKh4v8LyIRU9t5rQGGC4sCVlnTlQl4NWCqDrHx4_-iRi5L3AKVnjc7cSKwprfKPG4F994tXh2R1uWe66CwF1apt9PXa0hB2Ezz2_t_2Z55KtXkmwWNPGOCAJZDFKXt?key=AgyQUxUvRaYsce19bm41Iw)

Figure 8: FEA results of the upregulated genes of mutant glioma on LGr4

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcH0dGhdp94NF0hnfzfh_cUeyIhtX51adTVdWiIScu3t7KbuRuKBnW9lm_4NON_ITEctuba_4Ytys_yTkV0NaBP3ARiTeNYOLsw-gsM6WILxxtUfGGF2zVmmHOnEZDdOqFt5u1XAMNzTuhISj6O04s_TJk?key=AgyQUxUvRaYsce19bm41Iw)

Figure 9: FEA results of the downregulated genes of WT glioma on LGr1

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfue3EIIBdcOKeSBIsXrZDdQaQiIpqCpaVUr6h4Q4NuLFkpQw2P96G1rQCz-2nthDCIq7ViH7wrnNh_-Iy23qX6_0MfCn0ijCAiOwylkDtnl3vRPfo49bPSRtrE9yPv_1-Sw_qsAu7MkICOIZiwo-QnmBFF?key=AgyQUxUvRaYsce19bm41Iw)****

Figure 10: FEA results of the downregulated genes of WT glioma on LGr4

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcT_xyuhyU8OiOcH3rUiOfh_cRlF5noyBZnH79MC5PJSyGd6-52dmde9GsrMZj5SyIGOP7Ul3ARUE9bw6GPnjm8hZb3jxB8D5D9la63KW3fOijTaYO9MH_HFmdm-xMWzwVXNfvmC107LgiIUC0QGZq5G53y?key=AgyQUxUvRaYsce19bm41Iw)

Figure 11: FEA results of the downregulated genes of mutant glioma on LGr1

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdcsb_NghQE0mMNMeiOBmhbfL5K8jlhlmOV9dUtGFSywni5h2fw3_yitrNRN6RfHR4wCGWbekhQYPcIdcsxeDcjDRmODcpAgfs18CTyx2gdJDBYmkcuPKBDGCoATp8nhH9eHr0Qz9e7EJFj66r9r4PdL9ci?key=AgyQUxUvRaYsce19bm41Iw)

Figure 12: FEA results of the downregulated genes of mutant glioma on LGr4

**Comparison of results**

Our analysis showed that the up and down regulated pathways align with the pathways described in the paper. For example, some of the down regulated pathways in wildtype gliomas were related to cell adhesion, collagen metabolism, and extracellular matrix remodeling, which is consistent

with the poor prognosis and invasive behavior described in the paper. 

Additionally, immune-related pathways are significantly upregulated in mutant gliomas. Both

LGr1 and LGr4 analyses show enhanced lymphocyte activation, T-cell activation, and antigen

presentation via MHC class II, suggesting active immune responses

While the paper primarily focuses on the less aggressive nature of mutant gliomas, these immunoupregulations hint at complex interactions that affect immune response or tumor suppression.

**Machine Learning**

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfVwlsE3MkvY5Qr7sLaRlwIvS1p9UUyueE5qsb6HVNqnYD55VzeDid_LUm5GnYcfyHhK5XPWtMxlqm3baUdM__9T7Dx56I-S6ME2oERvCJ8l5-js3B1yAMfN6dHR062UMC2T93GPEaLKVEJLBpoBqa0Ui8l?key=AgyQUxUvRaYsce19bm41Iw)****

Figure 13: k-NN

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXd9Hc-PkjGcStu3NGVNcnZid58KaThQLT6_9GZbW1hOpRggUalA5HaEvu4GidL2Fq1Ju2LVkP3vZsI8JYMjbSIUDX-oWP9S0qcwoT7gJhJVvWRiEhT_wY-XNaYHqwz8snLAXJD6xZvuANRSKSZQg7E6mcE?key=AgyQUxUvRaYsce19bm41Iw)

Figure 14: Random Forest

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfDKAwtjDR_Ko2zbdS-RZtproYCLXeDAja05LyknYp63g-SCG-9uM_8XJmJLlXaAE3yG2lrFmIcIM7f_SsEHyegUaZMrfz2GmeyZY4oVGxH0mryUD7eWpGI1EQ941npF3xkOIUasAAD8BW-jpleSbYSWYbx?key=AgyQUxUvRaYsce19bm41Iw)

Figure 15: rf model



**Conclusion**

In summary, IDH-wildtype gliomas exhibit down regulated pathways linked to invasion and metastasis, reflecting their poor prognosis, while mutant gliomas show downregulation of proliferative pathways and upregulation of immune-related processes, consistent with slower growth and better clinical outcomes.

**Citations:**

1. Molinaro AM, Taylor JW, Wiencke JK, Wrensch MR. Genetic and molecular epidemiology of adult diffuse glioma. Nat Rev Neurol. 2019 Jul;15(7):405-417. doi: 10.1038/s41582-019-0220-2. Epub 2019 Jun 21. PMID: 31227792; PMCID:PMC7286557.

2. Whitfield BT, Huse JT. Classification of adult-type diffuse gliomas: Impact of the World Health Organization 2021 update. Brain Pathol. 2022 Jul;32(4):e13062. doi: 10.1111/bpa.13062. Epub 2022 Mar 14. PMID: 35289001;PMCID: PMC9245936.

3. Cairns R, Mak T. Oncogenic isocitrate dehydrogenase mutations: mechanisms, models, and clinical opportunities. Cancer Discovery 2013; 3: 730-41.

4. Ceccarelli et. al, (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell, 164(3), 550–563. <https://doi.org/10.1016/j.cell.2015.12.028>

**Word Count (Excluding in-text citations):** 419

