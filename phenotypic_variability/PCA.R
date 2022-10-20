load("rdata/phenotypes_regmap.RData")
load("rdata/phenotypes_toua.RData")
load("rdata/phenotypes_lang.RData")

library(ade4)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)
theme_set(theme_pubr())

regmap <- regmap[!str_detect(colnames(regmap), "ecotype")]
regmap$pop <- "regmap"
lang$pop <- "lang"
toua$pop <- "toua"

x <- rbind.data.frame(regmap, toua, lang)
x <- na.omit(x)
x <- x[,str_detect(colnames(x), "change|pop")]

pca <- dudi.pca(x[,!str_detect(colnames(x), "pop")], 
                scannf = FALSE, nf = 4, 
                center = TRUE, scale = TRUE)

pca$li$pop <- x$pop
pca$li$pop <- toupper(str_replace(pca$li$pop, "lang", "Languedoc"))
ade4::s.corcircle(pca$co, xax = 1, yax = 2) 
pca$co$name <- rownames(pca$co)

scree <- data.frame(component = seq(1:length(pca$eig)), 
                    eigen.values = pca$eig, 
                    explained.variance = round(pca$eig/sum(pca$eig) * 100, 2))

vars <- ggplot(pca$co, aes(x = Comp1, y = Comp2, 
                           label = str_split_fixed(rownames(pca$co), '_', 2)[,1], 
                           color = rownames(pca$co))) + 
  geom_segment(aes(xend=Comp1, yend=Comp2, x = 0, y = 0, 
                   color = rownames(pca$co)) , size = 1)+
  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  geom_point(size = 5)+ geom_label_repel(size=4)+
  scale_color_brewer(palette = "Set2") +
  xlab(paste("Principal component 1 :", 
             scree[1, "explained.variance"], "%"))+
  ylab(paste("Principal component 2 :", 
             scree[2, "explained.variance"], "%"))+
  theme(legend.position = "none") + 
  ggtitle("")

indivs <-ggplot(pca$li, aes(x = Axis1, y = Axis2, 
                            label = rownames(pca$li), color = pop)) + 
  geom_point(size = 2)+ stat_ellipse()+
  scale_color_brewer(palette = "Accent", name = "Population") + 
  ggtitle("") +
  xlab(paste("Principal component 1 :",
             scree[1, "explained.variance"], "%"))+
  ylab(paste("Principal component 2 :", 
             scree[2, "explained.variance"], "%"))

ggexport(indivs + vars+plot_annotation(tag_levels = 'a',), 
         filename = "results/pca_all_pop.pdf", width = 11)
