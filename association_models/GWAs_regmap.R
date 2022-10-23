# loads phenotype data
load("rdata/phenotypes_regmap.RData")

# loads genotype data (snp)
load("rdata/inferred_SNP_matrix_regmap.RData")

# loads snp locations in the genome
load("rdata/inferred_snp_chromosomic_map.RData")

library(statgenGWAS)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ade4)
library(stringr)

# individuals for which we have snp data
acc <- intersect(regmap$ecotype, colnames(snp))


# chromosomic positions
colnames(map) <- c("chr", "pos")
head(map)


# adding N/C ratio 
regmap$NC_aCo2 <- regmap$N_aCo2/regmap$C_aCo2
regmap$NC_eCo2 <- regmap$N_eCo2/regmap$C_eCo2
regmap$NC_change <- (regmap$NC_eCo2-regmap$NC_aCo2)/regmap$NC_aCo2*100


# adding PC1 
pca_change <- dudi.pca(na.omit(regmap[,str_detect(colnames(regmap), "change")]), 
                       scannf = FALSE, nf = 4, center = TRUE, scale = TRUE)
pca_aCo2 <- dudi.pca(na.omit(regmap[,str_detect(colnames(regmap), "aCo2")]), 
                     scannf = FALSE, nf = 4, center = TRUE, scale = TRUE)
pca_eCo2 <- dudi.pca(na.omit(regmap[,str_detect(colnames(regmap), "eCo2")]), 
                     scannf = FALSE, nf = 4, center = TRUE, scale = TRUE)
regmap$PC1_change <- pca_change$li[match(rownames(regmap), rownames(pca_change$li)),"Axis1"]
regmap$PC1_aCo2 <- pca_aCo2$li[match(rownames(regmap), rownames(pca_aCo2$li)),"Axis1"]
regmap$PC1_eCo2 <- pca_eCo2$li[match(rownames(regmap), rownames(pca_eCo2$li)),"Axis1"]


# formatting phenotypes
Y <- regmap[match(acc, regmap$ecotype), 
            stringr::str_detect(colnames(regmap), "change|eCo2|aCo2")]
Y$genotype <- acc
Y <- Y[,c("genotype", colnames(Y)[colnames(Y) != "genotype"])]


# formatting snps
X <- t(snp)


# kinship matrix : computes genetic relatedness between individuals
kin <- kinship(X, method = c("astle"))
#heatmap(kin)


# object storing all info for gwas analysis with statgenGWAS
gData <- createGData(geno = X, map = map, pheno = Y, kin = kin)


# removing low MAF snps
gData <- codeMarkers(gData, impute = FALSE, verbose = TRUE, MAF = 0.02) 
# Input contains 2678763 SNPs for 186 genotypes.
# 0 genotypes removed because proportion of missing values larger than or equal to 1.
# 0 SNPs removed because proportion of missing values larger than or equal to 1.
# 595198 SNPs removed because MAF smaller than 0.02.
# 1263152 duplicate SNPs removed.
# Output contains 820413 SNPs for 186 genotypes.
#summary(gData)

gData <- codeMarkers(gData, impute = FALSE, verbose = TRUE, MAF = 0.04)
# Input contains 820413 SNPs for 186 genotypes.
# 0 genotypes removed because proportion of missing values larger than or equal to 1.
# 0 SNPs removed because proportion of missing values larger than or equal to 1.
# 187719 SNPs removed because MAF smaller than 0.04.
# 0 duplicate SNPs removed.
# Output contains 632694 SNPs for 186 genotypes.


# to run a GWAS on a single phenotype, uncomment these lines :

# running mixed models for all snps and getting their pvalues
# gwas <- runSingleTraitGwas(gData = gData,
#                            traits = c("Cu_change"),
#                            GLSMethod = "single",
#                            remlAlgo = "EMMA",
#                            thrType = "fixed",
#                            LODThr = 6,
#                            nCores = 32)

# run for all elements at the same time (run independently, results merged)
gwas <- runSingleTraitGwas(gData = gData,
                           GLSMethod = "single",
                           remlAlgo = "EMMA",
                           nCores = 52, 
                           thrType = "small",
                           nSnpLOD = 50)



save(gwas, file = "rdata/gwas_all_traits_emma_top50_astle_maf0.04_mad5.rdata")



########## qqplots

qqs <- list()
for(trait in unique(t$trait)){
  qqs[[length(qqs)+1]] <- plot(gwas, plotType = "qq", trait = trait, output = F) + 
    ggtitle(trait)
}

gridExtra::grid.arrange(grobs = qqs)

ggexport(plotlist = qqs, 
         width = 1800, height = 1600, ncol = 5, nrow = 6,
         filename = "results/qqplots/all_qqs_maf0.04_regmap.png")


####################### Compare common SNPs between eCO2, aCO2 phenotypes
# and relative change

# 
# venns <- list()
# elements <- unique(stringr::str_split_fixed(colnames(regmap), '_', 2)[,1])
# for(el in elements[elements!= "ecotype"]){
#   DT <- as.data.table(gwas$signSnp$Y)
#   dt <- DT[stringr::str_detect(trait, paste0(el, '_'))]
#   venns[[el]] <- DIANE::draw_venn(list("relative change" = dt[dt$trait == paste0(el, '_', 'change'),]$snp,
#                                        "eCO2" = dt[dt$trait == paste0(el, '_', 'eCo2'),]$snp,
#                                        "aCO2" = dt[dt$trait == paste0(el, '_', 'aCo2'),]$snp)) + ggtitle(el)
# }
# gridExtra::grid.arrange(grobs = venns)

