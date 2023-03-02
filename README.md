# EC_metabolic_subtype_analysis
Metabolic subtype study for endomaterial cancer

## Workflow
<img src="https://github.com/Akankxha/EC_metabolic_subtype_analysis/blob/main/Figure1.jpg" width="600" />


## Available Scripts
`data_preprocessing.r` : Gene experssion data filtering and VST normalization \
`eda.r` : Exploratory data analysis including PCA \
`nmf_clustering.r` : NMF and Consensus clustering \
`clinical_profiles.r` : Survival analysis and fisher exact test \
`degs.r` : Differencial experssion analysis \
`prognostic_genes.r` : Univariate cox regression analysis \
`validation.r` : Validation of subtypes with GEO Datasets 

**Notes - ** This code repository contains work that are under review for a publication. Below are the details of the study.


### Identification and Characterization of Metabolic Subtypes of Endometrial Cancer using a Systems-Level Approach 
Endometrial cancer (EC) is the most common gynecological cancer worldwide. Understanding the metabolic adaptation and its heterogeneity in tumor tissues may provide new insights and help in cancer diagnosis, prognosis, and treatment. In this study, we investigated metabolic alterations of EC to understand the variations in the metabolism within tumor samples. Integration of tran-scriptomics data of EC (RNA-Seq) and the human genome-scale metabolic network was per-formed to identify the metabolic subtypes of EC and uncover the underlying dysregulated met-abolic pathways and reporter metabolites in each subtype. The relationship between metabolic subtypes and clinical variables was explored. Further, we correlated the metabolic changes oc-curring at the transcriptome level with the genomic alterations. Based on metabolic profile, EC patients are stratified into two subtypes (metabolic subtype-1 and subtype-2) that significantly correlate to patient survival, tumor stages, mutation, and copy number variations. We observed the co-activation of the pentose phosphate pathway, one-carbon metabolism, and genes involved in controlling estrogen levels in metabolic subtype-2, which is linked to poor survival. PNMT and ERBB2 are also upregulated in metabolic subtype-2 samples and present on the same chromo-some locus 17q12, which is amplified. PTEN and TP53 mutations show mutually exclusive be-havior between subtypes and display a difference in survival. This work identifies metabolic subtypes with distinct characteristics at the transcriptome and genome levels, highlighting the metabolic heterogeneity within EC.

**Keywords - ** metabolic reprogramming; endometrial cancer; transcriptome; systems biology; re-porter metabolites


**Question? ** For any query related to the code or the work, Please reach out to us at [akansha.srivastava@research.iiit.ac.in](mailto:akansha.srivastava@research.iiit.ac.in)
 