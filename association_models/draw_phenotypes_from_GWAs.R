library(ggpubr)
library(ggplot)

# loads the gwas object
load(file = "rdata/gwas_all_traits_emma_top50_astle_maf0.04_mad5.rdata")

# get significant SNPs
t <- gwas$signSnp$Y[,c("pValue", "chr", "pos", "allFreq", "snp", "trait")]
t <- t[order(t$pValue),]

d1 <- regmap[regmap$ecotype %in% rownames(X),]
d1[,c(t$snp)] <- 0
for(snip in t$snp){
  d1[,paste0( "snp", snip)] <- ifelse(X[as.character(d1$ecotype), snip] > 0, "1","0")
}

for(trait in unique(t$trait)){
  cond <- t$trait == trait
  snips <- paste0( "snp", t[cond,]$snp )                              
  d <- reshape2::melt(d1)
  d <- d[d$variable == trait,]
  plots <- list()
  for(snip in snips){
    plots[[length(plots)+1]] <- ggplot(d, aes_string(x = snip, y = "value", color = snip)) + 
      ggtitle(snip) +
      geom_jitter(size =1.2, width = 0.15, alpha = 0.6, color = "grey") + 
      scale_color_brewer(palette = "Set2") + 
      scale_fill_brewer(palette = "Set2") + 
      ggtitle(snip) + labs(subtitle = trait) + 
      geom_hline(yintercept = 0, size = 1, col = "grey") + 
      stat_summary(fun.data=median_q1q3, geom="pointrange", size = 0.9) + 
      theme_pubr(legend = "none") + labs_pubr()
  }
  ggexport(plotlist = plots, ncol = 10, width = 16000, height = 12000, nrow = 5,res = 600,
           filename = paste0("results/phenotype_", trait, "_top50_emma_astle_maf0.04.png"))
}
