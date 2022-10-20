# GWAS on the ionome response to elevated CO2

This repository contains code and partial data for a GWAs on mineral content variation under elevated CO2 in *Arabidopsis thaliana*.

Existing literature discloses that there is intra-specific variability in the phenotypic response to rising CO2 in Arabidopsis as well as in plants of agronomic interest. The variability observed in the ionome response to eCO2 was however never explained by genomic determinants in plants, even though it has the potential to fuel the discovery of new candidate genes controlling their mineral depletion. We detail here a GWAs project with the aim to identify genetic determinants of mineral status response to CO2 elevation. 


## Data

The `data` folder contains raw phenotype information about the three populations of Arabidopsis (Regmap, TOU-A, and Languedoc).


The `rdata` folder contains raw, pre-processed, and outliers-free phenotype phentoypes for the three populations.


## Pre-processing

The `data_preprocessing` folder contains several scripts for data import and cleaning:

+ `phenotypes_preprocessing.R` imports the raw csv phenotypes of all populations, computes the relative changes caused by eCO2 for all elements, and stores the result in the `rdata` folder.

+ `phenotypes_outlier_removal.R` imports the preprocessed phenotypes of all populations, and replaces observations that are more than 5 median absolute deviations away from the median by NAs for each phenotype in each population. It  stores the results in the `rdata` folder.



## Analysis of the natural variability in three populations of Arabidopsis


The `phenotypic_variability` folder contains the analyses of the natural variability in three populations of Arabidopsis.


+ `phenotypes_distributions.R` shows the distribution of all relative changes in the three populations and exports the plot in the `results`.

+ `PCA.R` shows the PCA analysis of the relative changes in the three populations and exports the plot in the `results`.

+ `clustering_regmap.R` applies a k-means clustering to the accessions of the Regmap panel based on their relative changes and identifies a group of resilisent accessions. It exports the clustering plot in the `results`.