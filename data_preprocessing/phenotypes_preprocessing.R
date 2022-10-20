library(stringr)


#################################### REGMAP ######################################

# data import and saving to rdata

regmap <- read.csv("data/phenotypes_regmap.csv", na.strings = ".")

# Hov 4 est dupliqué, un des seuls réplicats peut être. J'en fais une moyenne
hov4.1 <- colMeans(regmap[regmap$genotype == "Hov4-1", colnames(regmap) != "genotype"])

regmap <- regmap[regmap$genotype != "Hov4-1",]
rownames(regmap) <- regmap$genotype
regmap <- regmap[,colnames(regmap) != "genotype"]
regmap["Hov4-1",] <- hov4.1
colnames(regmap)

annot <- read.csv("data/accessions_1307_regMap.csv")
#save(annot, file = "rdata/accessions_regmap_annotation.RData")

rownames(regmap)[!rownames(regmap) %in% annot$original_n..]

# ids from 1001 genome project
regmap$ecotype <- annot[match(rownames(regmap), annot$original_n..), "Ecotype_ID"]

elements <- unique(stringr::str_split_fixed(colnames(regmap), '_', 2)[,1])
elements <- elements[elements != "ecotype"]

for(el in elements){
  regmap[,paste0(el, '_change')] <- (regmap[,paste0(el,"_eCo2")] - regmap[,paste0(el,"_aCo2")])/
    regmap[,paste0(el,"_aCo2")] * 100
}


save(regmap, file = "rdata/phenotypes_regmap_unfiltered.RData")


#################################### LANGUEDOC ######################################

# data import and saving to rdata


lang <- read.csv("data/LANGUEDOC_phenotypes.csv", na.strings = ".", sep = ';')
rownames(lang) <- lang$genotype


lang <- lang[,colnames(lang) != "genotype"]
lang <- lang[,!str_detect(colnames(lang), "ratio|CN")]

elements <- unique(stringr::str_split_fixed(colnames(lang), '_', 2)[,1])
elements <- elements[elements != "ecotype"]

for(el in elements){
  lang[,paste0(el, '_change')] <- (lang[,paste0(el,"_eCo2")] - lang[,paste0(el,"_aCo2")])/
    lang[,paste0(el,"_aCo2")] * 100
}

save(lang, file = "rdata/phenotypes_lang_unfiltered.RData")




#################################### TOUA ######################################

# data import and saving to rdata

toua <- read.csv("data/TOUA_phenotypes.csv", na.strings = ".", sep = ';')
rownames(toua) <- toua$genotype


toua <- toua[,colnames(toua) != "genotype"]
toua <- toua[,!str_detect(colnames(toua), "ratio|CN")]

elements <- unique(stringr::str_split_fixed(colnames(toua), '_', 2)[,1])

for(el in elements){
  toua[,paste0(el, '_change')] <- (toua[,paste0(el,"_eCo2")] - toua[,paste0(el,"_aCo2")])/
    toua[,paste0(el,"_aCo2")] * 100
}

save(toua, file = "rdata/phenotypes_toua_unfiltered.RData")
