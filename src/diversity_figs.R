library(tidyverse)
library(plotly)
library(htmlwidgets)
library(colorspace)
setwd("~/WorkForaging/Academia/Nicole/")

###Alpha Diversity

#shannon alpha diversity
data <- read_tsv('data/diversity/exported/shannon_vector/alpha-diversity.tsv') %>% rename(SRAid= '...1')
meta <- read_csv('data/SraRunTable.csv') %>% rename(SRAid = "Run")

head(data)
gdat <- left_join(meta, data, by='SRAid') %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

ggplot(gdat)+
  geom_boxplot(mapping=aes(x=Shrub, y=shannon_entropy, fill=Shrub))+
  geom_jitter(mapping=aes(x=Shrub, y=shannon_entropy), width = 0.2)+
  ylab("Shannon Entropy")+
  facet_wrap(~Site)+
  theme_bw()+
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, vjust = .7))
ggsave('figs/shannon_entropy.png', width=4000, height=3000, units='px', dpi=300)

#peilou's evenness (shannon / # species)
data <- read_tsv('data/diversity/exported/evenness_vector/alpha-diversity.tsv')%>% rename(SRAid= '...1')
meta <- read_csv('data/SraRunTable.csv') %>% rename(SRAid = "Run")

gdat <- left_join(meta, data, by='SRAid') %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

ggplot(gdat)+
  geom_boxplot(mapping=aes(x=Shrub, y=pielou_evenness, fill=Shrub))+
  geom_jitter(mapping=aes(x=Shrub, y=pielou_evenness), width = 0.2)+
  ylab("Pielou Evenness")+
  facet_wrap(~Site)+
  theme_bw()+
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, vjust = .6))

ggsave('figs/pielou_evenness.png', width=4000, height=3000, units='px', dpi=300)

#faith's phylogenetic diversity
data <- read_tsv('data/diversity/exported/faith_pd_vector/alpha-diversity.tsv')%>% rename(SRAid= '...1')
meta <- read_csv('data/SraRunTable.csv') %>% rename(SRAid = "Run")

gdat <- left_join(meta, data, by='SRAid') %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

ggplot(gdat)+
  geom_boxplot(mapping=aes(x=Shrub, y=faith_pd, fill=Shrub))+
  geom_jitter(mapping=aes(x=Shrub, y=faith_pd), width = 0.2)+
  ylab("Faith's Phylogenetic Diversity")+
  facet_wrap(~Site)+
  theme_bw()+
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, vjust = .6))
ggsave('figs/faith_pd.png', width=4000, height=3000, units='px', dpi=300)


####Beta Diversity

#bray curtis distance matrix
data <- read_tsv('data/diversity/exported/bray_curtis_pcoa_results/ordination.txt',
                 skip=9, n_max=50, col_names = c('SRAid', paste("PC", 1:51, sep='')))

gdat <- left_join(meta, data %>% select(SRAid, PC1, PC2, PC3), by="SRAid") %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

plt <- plot_ly(data=gdat, x=~PC1, y=~PC2, z=~PC3, color=~Shrub, symbol=~Site, 
        symbols = c("circle", "square", "diamond"), type = 'scatter3d',mode='markers') %>%
  layout(title = "Bray-Curtis Distance PCoA",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ))
plt
saveWidget(plt, 'figs/bray_curtis_pca.html', selfcontained = T)

#jaccard distance
data <- read_tsv('data/diversity/exported/jaccard_pcoa_results/ordination.txt',
                 skip=9, n_max=50, col_names = c('SRAid', paste("PC", 1:51, sep='')))
gdat <- left_join(meta, data %>% select(SRAid, PC1, PC2, PC3), by="SRAid") %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

plt <- plot_ly(data=gdat, x=~PC1, y=~PC2, z=~PC3, color=~Shrub, symbol=~Site, 
               symbols = c("circle", "square", "diamond"), type = 'scatter3d',mode='markers') %>%
  layout(title = "Jaccard Distance PCoA",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")
         ))
plt
saveWidget(plt, 'figs/jaccard_pca.html', selfcontained = T)

#unweighted_unifrac distance
data <- read_tsv('data/diversity/exported/unweighted_unifrac_pcoa_results/ordination.txt',
                 skip=9, n_max=50, col_names = c('SRAid', paste("PC", 1:51, sep='')))
gdat <- left_join(meta, data %>% select(SRAid, PC1, PC2, PC3), by="SRAid") %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

plt <- plot_ly(data=gdat, x=~PC1, y=~PC2, z=~PC3, color=~Shrub, symbol=~Site, 
               symbols = c("circle", "square", "diamond"), type = 'scatter3d',mode='markers') %>%
  layout(title = "Unweighted Unifrac Distance PCoA",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")
         ))
plt
saveWidget(plt, 'figs/unweighted_unifrac_pca.html', selfcontained = T)

#weighted_unifrac distance
data <- read_tsv('data/diversity/exported/weighted_unifrac_pcoa_results/ordination.txt',
                 skip=9, n_max=50, col_names = c('SRAid', paste("PC", 1:51, sep='')))
gdat <- left_join(meta, data %>% select(SRAid, PC1, PC2, PC3), by="SRAid") %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush")))

plt <- plot_ly(data=gdat, x=~PC1, y=~PC2, z=~PC3, color=~Shrub, symbol=~Site, 
               symbols = c("circle", "square", "diamond"), type = 'scatter3d',mode='markers') %>%
  layout(title = "Weighted Unifrac Distance PCoA",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")
         ))
plt
saveWidget(plt, 'figs/weighted_unifrac_pca.html', selfcontained = T)
