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


####################### Starting some analyses for GWAs quality

# get significant SNPs
t <- gwas$signSnp$Y[,c("pValue", "chr", "pos", "allFreq", "snp", "trait")]
t <- t[order(t$pValue),]

########## qqplots

qqs <- list()
for(trait in unique(t$trait)){
  qqs[[length(qqs)+1]] <- plot(gwas, plotType = "qq", trait = trait, output = F) + 
    ggtitle(trait)
}

gridExtra::grid.arrange(grobs = qqs)

ggexport(plotlist = qqs, 
         width = 1800, height = 1600, ncol = 5, nrow = 6,
         filename = "results/all_qqs_maf0.04_regmap.png")


###########  Variance plots

d <- data.frame(gwas$GWASInfo$varComp$Y)
d$var <- rownames(d)
d <- reshape2::melt(d)
ggplot(d, aes(x = variable, y = value, 
              label = round(value, 1), 
              fill = var)) + geom_bar(position="stack", stat = "identity") + 
  facet_wrap(~variable, scales = "free") + 
  ggtitle("Variance estimates for mixed models random effects") + 
  labs(subtitle = "Vg = genetic variance, Ve = residual variance") +
  scale_fill_brewer(palette = "Set2") + theme_pubr()


############# allelic frequency of snps :

ggplot(t[order(t$allFreq),], 
       aes(x = snp, size = allFreq, color = trait, 
           y = round(allFreq, 4))) + 
  geom_point() + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + coord_flip()



