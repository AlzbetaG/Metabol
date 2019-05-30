# Metabol
The package of preprocessing and statistical analysis of metabolomic data sets made by logratio methodology.

[![DOI](https://zenodo.org/badge/188906403.svg)](https://zenodo.org/badge/latestdoi/188906403)

What is it?
-----------

- The preprocessing: quality control samples treatment and zero imputation
- The comparison of the logratio methodology with the other transformations/scalings: log, ln, PQN, no transformation, ln(PQN), Pareto 
scaling
- Multivariate methods: PCA, PLS-DA (also with VIP plots and permutation test), VIP plots, OPLS-DA (also with permutation test), cluster analysis
- Univariate methods: boxplots, ROC curves, volcano plots
- Statistical tests: Shapiro-Wilk test, ANOVA, Kruskal-Wallis test, t-test, permutation test 
- The analysis of correlations in data
- The analysis of outliers in data

Getting Started
---------------

### Dependencies

The package has dependencies on 

	R (>= 2.10), openxlsx, sm, ropls

### Installation

Installion of `Metabol` is really easy for registered users (when the R-tools are installed). Just use 

    library(devtools)
    install_github("Metabol", "AlzbetaG")
