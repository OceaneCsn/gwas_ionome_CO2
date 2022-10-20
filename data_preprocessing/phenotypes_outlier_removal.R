load("rdata/phenotypes_regmap_unfiltered.RData")
THR <- 5

primary_data <- regmap[!stringr::str_detect(colnames(regmap), "change|ecotype")] 

for(col in colnames(primary_data)){
  
  vals <- primary_data[,col]
  med <- median(na.omit(vals))
  sd <- mad(na.omit(vals))
  
  condition <- vals < med - THR*sd | vals > med + THR*sd 
  primary_data[,col] <- ifelse(condition, NA, primary_data[,col])
}

elements <- unique(stringr::str_split_fixed(colnames(primary_data), '_', 2)[,1])

for(el in elements){
  primary_data[,paste0(el, '_change')] <- (primary_data[,paste0(el,"_eCo2")] - primary_data[,paste0(el,"_aCo2")])/
    primary_data[,paste0(el,"_aCo2")] * 100
}

primary_data$ecotype <- regmap$ecotype 
regmap <- primary_data


save(regmap, file = "rdata/phenotypes_regmap.RData")


######################## LANGUEDOC

load("rdata/phenotypes_lang_unfiltered.RData")
primary_data <- lang[!stringr::str_detect(colnames(lang), "change|ecotype")] 
for(col in colnames(primary_data)){
  
  vals <- primary_data[,col]
  med <- median(na.omit(vals))
  sd <- mad(na.omit(vals))
  
  condition <- vals < med - THR*sd | vals > med + THR*sd 
  primary_data[,col] <- ifelse(condition, NA, primary_data[,col])
}

elements <- unique(stringr::str_split_fixed(colnames(primary_data), '_', 2)[,1])

for(el in elements){
  primary_data[,paste0(el, '_change')] <- (primary_data[,paste0(el,"_eCo2")] - primary_data[,paste0(el,"_aCo2")])/
    primary_data[,paste0(el,"_aCo2")] * 100
}
lang <- primary_data
save(lang, file = "rdata/phenotypes_lang.RData")



######################## TOUA

load("rdata/phenotypes_toua_unfiltered.RData")
primary_data <- toua[!stringr::str_detect(colnames(toua), "change|ecotype")] 
for(col in colnames(primary_data)){
  
  vals <- primary_data[,col]
  med <- median(na.omit(vals))
  sd <- mad(na.omit(vals))
  
  condition <- vals < med - THR*sd | vals > med + THR*sd 
  primary_data[,col] <- ifelse(condition, NA, primary_data[,col])
}

elements <- unique(stringr::str_split_fixed(colnames(primary_data), '_', 2)[,1])

for(el in elements){
  primary_data[,paste0(el, '_change')] <- (primary_data[,paste0(el,"_eCo2")] - primary_data[,paste0(el,"_aCo2")])/
    primary_data[,paste0(el,"_aCo2")] * 100
}
toua <- primary_data
save(toua, file = "rdata/phenotypes_toua.RData")
