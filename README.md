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

+ `regmap_imputed_SNPs_preprocessing.R` imports the raw imputed SNP data from Arouisse et Al. (https://figshare.com/articles/dataset/arabidopsis_2029_Maf001_filter95/11346875), that should be placed in the `data` folder (data not included to this repository because of their size). It reads it with appropriate tools (**This step is requires a lot of memory (RAM)**), keeps only the accessions of the regmap panel present in our experiment, and removes SNPs that do not occurr in the accessions of our experiment. It stores the results in the `rdata` folder.


## Analysis of the natural variability in three populations of Arabidopsis


The `phenotypic_variability` folder contains the analyses of the natural variability in three populations of Arabidopsis.


+ `phenotypes_distributions.R` shows the distribution of all relative changes in the three populations and exports the plot in the `results`.

+ `PCA.R` shows the PCA analysis of the relative changes in the three populations and exports the plot in the `results`.

+ `clustering_regmap.R` applies a k-means clustering to the accessions of the Regmap panel based on their relative changes and identifies a group of resilisent accessions. It exports the clustering plot in the `results`.


## Association models on the REGMAP panel


+ `GWAs_regmap.R` runs Linear Mixed Models on the Regmap Panel via the StatgenGWAS R package. It first cleans the SNP matrix by removing low MAF SNPs and duplicated SNPs, and then runs the association models on each phenotype. Stores the results in the `rdata` folder. Also draws quantile-quantile (qq) plots for GWAs quality control, as well as variance decomposition and allelic frequencies of top 50 SNPs.

+ `draw_manhattan_from_GWAs.R` plots the Manhattan plots of the GWAs for each phenotypes. 

+ `export_manhattan_from_GWAs_for_IGV.R` writes the Manhattan plots of the GWAs for each phenotypes, so that they can directly be imported in the IGV software and interactively browsed along the Arabidopsis genome.

+ `draw_phenotypes_from_GWAs.R` plots the phenotype of plants depending on weather they possess the top 50 SNPs of the GWAs or not, thus showing the effect of top variants on phenotypic responses. 
