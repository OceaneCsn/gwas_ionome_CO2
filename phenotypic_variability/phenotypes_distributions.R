library(ggpubr)
library(patchwork)
library(ggplot2)
library(stringr)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R") 


load("rdata/phenotypes_regmap.RData")
d <- reshape2::melt(regmap[,stringr::str_detect(colnames(regmap), "change")])
d$variable <- str_remove(d$variable, '_change')

dist <- ggplot(d, aes(x = variable, y = value, fill = variable)) + 
  geom_hline(yintercept = 0, size = 2, col = "darkred") +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.9) + 
  geom_point(aes(y = value, color = variable), position = position_jitter(width = 0.15), size = 1.25, alpha = 0.9) + 
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, fill = "white") + 
  labs(y = "Relative change (%)", x = NULL) + 
  ggtitle("REGMAP") + 
  scale_color_brewer(palette = "Set2", name = "Element")+ 
  scale_fill_brewer(palette = "Set2", name = "Element") + 
  theme(axis.text = element_text(size=16), axis.title = element_text(size=15))+ 
  theme_pubr(legend = "none")+ ylim(-70, 150)


load("rdata/phenotypes_toua.RData")
d <- reshape2::melt(toua[,stringr::str_detect(colnames(toua), "change")])
d$variable <- str_remove(d$variable, '_change')


toua <- ggplot(d, aes(x = variable, y = value, fill = variable)) + 
  geom_hline(yintercept = 0, size = 2, col = "darkred") +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.9) + 
  geom_point(aes(y = value, color = variable), position = position_jitter(width = 0.15), size = 1.25, alpha = 0.9) + 
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, fill = "white") + 
  labs(y = "Relative change (%)", x = NULL) + 
  ggtitle("TOU-A") + 
  scale_color_brewer(palette = "Set2", name = "Element")+ 
  scale_fill_brewer(palette = "Set2", name = "Element") + 
  theme_pubr(legend = "none")+ ylim(-70, 150)


load("rdata/phenotypes_lang.RData")
d <- reshape2::melt(lang[,stringr::str_detect(colnames(lang), "change")])
d$variable <- str_remove(d$variable, '_change')


lang <- ggplot(d, aes(x = variable, y = value, fill = variable)) + 
  geom_hline(yintercept = 0, size = 2, col = "darkred") +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.9) + 
  geom_point(aes(y = value, color = variable), position = position_jitter(width = 0.15), size = 1.25, alpha = 0.9) + 
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, fill = "white") + 
  labs(y = "Relative change (%)", x = NULL) + 
  ggtitle("LANGUEDOC") + 
  scale_color_brewer(palette = "Set2", name = "Element")+ 
  scale_fill_brewer(palette = "Set2", name = "Element") + 
  theme_pubr(legend = "none")+ ylim(-70, 150)

fig <-  dist + lang + toua + plot_annotation(tag_levels = 'a') ; fig

ggexport(fig,
         filename = "results/phenotype_distributions.pdf",
         width = 13, height = 5)

