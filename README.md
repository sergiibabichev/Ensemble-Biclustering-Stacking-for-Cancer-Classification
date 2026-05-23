# \# Biclustering-Driven Feature Engineering for Robust Cross-Platform Generalization in Transcriptomic Cancer Classification

# This repository contains the implementation of the ensemble biclustering and stacking framework proposed in the manuscript:

# \*\*“Biclustering-Driven Feature Engineering for Robust Cross-Platform Generalization in Transcriptomic Cancer Classification”\*\*

# The project focuses on robust feature engineering for high-dimensional transcriptomic cancer data using ensemble biclustering, consensus integration, and classifier-independent validation.

# \---

# \# Repository Structure

# \## Ensemble Biclustering Algorithms

# \### BCCC-based Ensemble

# \- `BCCC\_Ensemble.R`  

# &#x20; Implementation of the ensemble BCCC biclustering framework.

# \- `BCCC\_Ensemble\_Analysis.R`  

# &#x20; Structural and biological analysis of BCCC-derived biclusters.

# \- `BCCC\_Parameters optimization.R`  

# &#x20; Hyperparameter optimization for the BCCC algorithm.

# \---

# \### ISA-based Ensemble

# \- `ISA\_Ensemble.R`  

# &#x20; Implementation of the ensemble ISA biclustering framework.

# \- `ISA\_Ensemble\_Analysis.R`  

# &#x20; Structural and biological analysis of ISA-derived biclusters.

# \- `ISA\_Parameters optimization.R`  

# &#x20; Hyperparameter optimization for the ISA algorithm.

# \---

# \### PMD-based Ensemble

# \- `PMD\_Ensemble.R`  

# &#x20; Implementation of the ensemble PMD biclustering framework.

# \- `PMD\_Ensemble\_Analysis.R`  

# &#x20; Structural and biological analysis of PMD-derived biclusters.

# \- `PMD\_Parameters optimization.R`  

# &#x20; Hyperparameter optimization for the PMD algorithm.

# \---

# \## Stacking Integration

# \- `Stacking.R`  

# &#x20; Full stacking framework integrating BCCC + ISA + PMD consensus matrices.

# \- `Stacking\_short.R`  

# &#x20; Reduced stacking framework integrating ISA + PMD consensus matrices.

# \---

# \## Classification and Validation

# \### Random Forest and XGBoost Evaluation

# \- `cancer\_BCCC\_RF\_XGB\_Classification.py`

# \- `cancer\_ISA\_RF\_XGB\_Classification.py`

# \- `cancer\_PMD\_RF\_XGB\_Classification.py`

# \- `cancer\_STACK\_full\_RF\_XGB\_Classification.py`

# \- `cancer\_STACK\_reduced\_RF\_XGB\_Classification.py`

# These scripts perform:

# \- feature subset loading,

# \- train/test splitting,

# \- Bayesian hyperparameter optimization,

# \- Random Forest classification,

# \- XGBoost classification,

# \- cross-validation,

# \- independent validation,

# \- evaluation using Accuracy and Weighted F1-score.

# \---

# \# Methodological Overview

# The proposed framework includes:

# 1\. Preprocessing and filtering of transcriptomic gene expression data.

# 2\. Ensemble biclustering using:

# &#x20;  - ISA,

# &#x20;  - PMD,

# &#x20;  - BCCC.

# 3\. Consensus matrix construction from bootstrap realizations.

# 4\. Representation-level stacking of consensus masks.

# 5\. Greedy extraction of stable biclusters.

# 6\. Biological validation using GO and KEGG enrichment analysis.

# 7\. Classification using Random Forest and XGBoost classifiers.

# \---

# \# Main Objectives

# The framework is designed to:

# \- improve stability of gene selection,

# \- preserve biologically coherent co-expression structures,

# \- enhance cross-platform robustness,

# \- support interpretable transcriptomic cancer classification.

# \---

# \# Data

# The study was performed using:

# \- TCGA RNA-Seq dataset GEO GSE304485,

# \- GEO GSE2109 (expO) microarray dataset.

# Processed datasets are publicly available in GEO:

# \- \*\*GSE304485\*\*

# \---


# \# Technologies and Libraries

# \## R

# Main packages include:

# \- biclust

# \- isa2

# \- PMA

# \- clusterProfiler

# \- org.Hs.eg.db

# \- TCGAbiolinks

# \- GEOquery

# \- ggplot2

# \## Python

# Main libraries include:

# \- scikit-learn

# \- xgboost

# \- pandas

# \- numpy

# \- matplotlib

# \- bayesian-optimization

# \---

# \# Authors

# \- Sergii Babichev



