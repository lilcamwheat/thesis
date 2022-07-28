# MPhil Thesis

This repo contains all the code used for my MPhil Thesis. The data contains patient information and cannot be made publically available but this code can be used with other metabolomics. The primary package used in this project was the lilikoi R package (https://pubmed.ncbi.nlm.nih.gov/33484242/). An additional package used was the survey package (https://cran.r-project.org/web/packages/survey/survey.pdf), which was useful to adapt the prognosis model to case-cohort studies.

## Main Code

The commented code used for most of the analysis (data wrangling, pre-processing, exploratory analysis, testing proportional hazards, downstream pathway analysis, correlation heatmap) can be found under the main_code.R document.

## Edited lilikoi code to retrieve variable importance scores

The lilikoi package does not output importance scores for features for different methods. As such, I implemented changes to the source code package to extract the variable importance scores for machine learning. The code with the changes can be found under the edited_ML.R file. This code trains the classification models with 10CV and 100repeats, extracts the most important features and generates upset plots.

## Plotting

The code used to generate Figures 9 and 14, as well as Supplementary Figures 5 and 8 can be found under the plotting.R file.

## Prognosis prediction and Cross-validated KM curves function

The km_cv_ccs.R function fits a prognosis model, Cox-PH, to the data while accounting for the case-cohort design of the study. It also generates cross-validated KM curves (Figures 10, 11, and 15, as well as Supplemenrary Figures 3 and 9)

## Permutation test function

The permutation_ccs.R function performs 500 permutations to identify the significance level for the log-rank test of the cross-validated KM curves for case-cohort studies.
