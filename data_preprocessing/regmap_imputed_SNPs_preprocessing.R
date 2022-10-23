library("genio")

file_plink <- 'data/arabidopsis_2029_Maf001_filter95'

############ attention
## il faut 45Go de RAM pour faire cette opération :')
## fait sur le serveur plasticite
 
time_read_genio <- system.time(
  data_genio <- read_plink(file_plink)
)

############################################################# SNPs values


# X contient les 2029 accessions pour 3 millions de variants environ
# 2029 = les 894 uniques à la regmap + les 1135 de 1001 génomes

dim(data_genio$X)
x <- data_genio$X
head(colnames(x))
head(rownames(x))

# on commence par garder uniquement les accessions dont on va avoir besoin
load("rdata/phenotypes_regmap.RData")
common_acc <- intersect(as.character(regmap$ecotype), colnames(x))
# on retrouve bien nos 186 accessions dans la matrice de génotypes inférés

snp <- data.frame(x[,common_acc], check.names = F)
snp[1:10,1:10]
dim(snp)
snp <- snp[rowSums(snp) > 0,]
save(snp, file = "rdata/inferred_SNP_matrix_regmap.RData")

############################################################# SNPs locations


# snp locations in the genome
head(data_genio$bim)
map <- data.frame(data_genio$bim)
rownames(map) <- map$id
map <- map[,c("chr", "pos")]
colnames(map) <- c("chrom", "pos")
map <- map[rownames(snp),]
save(map, file = "rdata/inferred_snp_chromosomic_map.RData")
tail(rownames(snp))