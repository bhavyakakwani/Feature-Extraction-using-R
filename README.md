# Feature Extraction using R

Project by: **Bhavya Kakwani**

***

# Aim

Write a R program that extracts the features (aliphatic index, Boman index, hmoment index, peptide charge, etc.) of each “Peptide Sequence” given in the input file and create a feature matrix.

# Packages Used

1) Peptides

2) peptider

# Input Format

Through Command Line:

```r
rscript feature_extraction.R input.csv
```

# Ouput Format

Generates the following files upon execution of the R script:

1) **output.csv**: The file with the required feature matrix.

2) **logs.txt**: Any exceptions/errors are logged in this ffile.