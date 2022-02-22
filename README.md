# sandu-2015

An attempt to reproduce the mRNA and miRNA expression matrices and raw data preprocessing pipelines from the following publication: 

Sandhu, V., Bowitz Lothe, I. M., Labori, K. J., Lingjærde, O. C., Buanes, T., Dalsgaard, A. M., … Kure, E. H. (2015). Molecular signatures of mRNAs and miRNAs as prognostic biomarkers in pancreatobiliary and intestinal types of periampullary adenocarcinomas. Molecular Oncology, 9(4), 758–771. https://doi.org/10.1016/J.MOLONC.2014.12.002

The raw data are publicly available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60980 (GSE60980_RAW.tar)

# Reproducibility Notes

In general we couldn't get exactly the same matrices, but we did get some using diffent methods from the [AgiMicroRna](https://bioconductor.org/packages/release/bioc/html/AgiMicroRna.html) `R` package (v2.36).

Some notes:

- `reproduce.R`: use MRA normalization on both mRNA and miRNA data seems to be the best option in terms of number of produced features
- `select_var_features.R`: remove features that are not variable enough (as was done in the paper)

