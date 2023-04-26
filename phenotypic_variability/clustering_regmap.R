library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(patchwork)

load("rdata/phenotypes_regmap.RData")
# annot <- read.csv("data/accessions_1307_regMap.csv")

x <- na.omit(regmap)
x <- x[,str_detect(colnames(x), "change")]


set.seed(125)
kmeans <- kmeans(x, centers = 3)
x$accession <- rownames(x)
d <- reshape2::melt(x)
d$cluster <- kmeans$cluster[d$accession]

table(kmeans$cluster)

fig <- ggplot(d, aes(x = variable, y = value, fill = factor(cluster)))+ 
  geom_hline(yintercept = 0, size = 1, col = "darkred") + 
  geom_boxplot(width = 0.5) + 
  scale_fill_brewer(palette = "Accent", name = "Cluster") + 
  ggtitle("Cluster profiles in mineral status response")+ 
  ylab("Change (%)") + 
  scale_x_discrete(breaks = unique(d$variable),
                   labels = str_split_fixed(unique(d$variable), '_', 2)[,1])+
  xlab("")+theme_pubr()+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))+
  ggplot(d, aes(x = variable, y = value, fill = factor(cluster), color = factor(cluster))) + 
  geom_line(alpha = 0.5, size=1, aes(group = accession, color = factor(cluster))) +
  ggtitle("Cluster profiles in mineral status response") + 
  scale_color_brewer(name = "Cluster", palette = "Accent") +
  scale_fill_brewer(name = "Cluster", palette = "Accent") +
  geom_dotplot(color = 'black', binaxis='y', stackdir='center', 
               dotsize = 100, binwidth = 1/40)+
  geom_hline(yintercept = 0, size = 1, col = "darkred") + 
  scale_x_discrete(breaks = unique(d$variable),
                   labels = str_split_fixed(unique(d$variable), '_', 2)[,1])+
  xlab("") + ylab("Change (%)")+theme_pubr()+ 
  facet_wrap(~cluster, nrow = 1)+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  plot_layout(widths = c(0.7,1)) + plot_annotation(tag_levels = 'a');fig

ggexport(fig, filename = "results/clustering_regmap.pdf", width = 15, height = 7)


## cluster profiles only:

profiles <- ggplot(d, aes(x = variable, y = value, fill = factor(cluster), color = factor(cluster))) + 
  geom_line(alpha = 0.5, size=1, aes(group = accession, color = factor(cluster))) +
  ggtitle("Cluster profiles in mineral status response") + 
  scale_color_brewer(name = "Cluster", palette = "Accent") +
  scale_fill_brewer(name = "Cluster", palette = "Accent") +
  geom_dotplot(color = 'black', binaxis='y', stackdir='center', 
               dotsize = 100, binwidth = 1/40)+
  geom_hline(yintercept = 0, size = 1, col = "darkred") + 
  scale_x_discrete(breaks = unique(d$variable),
                   labels = str_split_fixed(unique(d$variable), '_', 2)[,1])+
  xlab("") + ylab("Change (%)")+theme_pubr()+ 
  facet_wrap(~cluster, nrow = 1)+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12)) 


ggexport(profiles, filename = "results/clustering_regmap_profiles_only.pdf", width = 10, height = 7)


## list of accessions in all clusters
acc_list <- d %>%
  select(accession, cluster) %>%
  group_by(accession) %>%
  reframe(cluster = cluster) %>%
  distinct() %>%
  arrange(cluster)

write.table(acc_list, quote = F, row.names = F, sep = "\t", 
            file = "results/accessions_per_clusters.tsv")
