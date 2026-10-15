library(phyloseq)
library(dplyr)
library(tidyverse)
library(vegan)
library(decontam)
library(picante)
setwd("~/WorkForaging/Academia/Nicole/drhizo/")

#import data
shrublevs <- c("Creosote","Mariola","Mesquite","Tarbush","Interspace")
shrubpal <- c("#ED9F56","#7AC2CF","#A479A5","#5EBA49","#EFD273")

ps <- readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds")
tax <- data.frame(tax_table(ps))
data <- data.frame(t(otu_table(ps))) 
meta <- data.frame( sample_data(ps) )  %>% rownames_to_column("SampleID")

#format for vegan
otu_id <- colnames(data)
asv <- otu_table(as.matrix(data), taxa_are_rows = F)
rownames(asv) <- rownames(data)
mdat <- sample_data(meta %>% 
                      select(SampleID, Site, Shrub) %>% 
                      mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush"))))
rownames(mdat) <- mdat$SampleID
ttable <- tax_table(as.matrix(tax %>% select(setdiff(colnames(tax), c("OTU","Taxon","Confidence")))) )
rownames(ttable) <- rownames(tax)
ps <- phyloseq(asv, mdat, ttable)

#filter data
conc <- colSums(asv)
contamdf <- isContaminant(ps, method="frequency",conc=conc, threshold=.1)
sum(contamdf$contaminant) 
ps_filt <- prune_taxa(!contamdf$contaminant, ps) #57 ASV removed at .1 thresh
ps_filt <- subset_samples(ps_filt, SampleID != 'SRR14874590' & SampleID != 'SRR14874591') #2 samples removed for low counts
ps_filt <- subset_taxa(ps_filt, Domain == "Bacteria") #filter for only Bacteria; 74 taxa removed 

#rarefy
otu <- as(otu_table(ps), "matrix")
otu <- as.data.frame(otu)
vegan::rarecurve(otu, step = 100, cex = 0.75, las = 1, label = F) #???
conc #see lib sizes
rare <- rrarefy(otu, sample = 10000) #seeing some disagreement online by the maintainers of vegan that this may be unhelpful
ps_r <- phyloseq(otu_table(rare, taxa_are_rows = F) ) #, subset_samples(sample_data( meta %>% rownames_to_column("SampleID"))  , SampleID %in% rownames(rare)))

#calculate alpha richness
richness <- estimate_richness(ps_r)
richness <- estimate_richness(ps_r, measures = c("Observed","Shannon","Simpson"))
anova.shrub.site = aov(richness$Shannon ~ sample_data(ps_filt)$Shrub + sample_data(ps_filt)$Site)
summary(anova.shrub.site)

#graph it up
gdat <- richness
gdat$SampleID <- rownames(gdat)
gdat <- left_join(gdat, data.frame(sample_data(ps)), by='SampleID') %>% 
  mutate( Evenness = Shannon / log(Observed)) %>% 
  pivot_longer(all_of(c('Observed','Shannon','Simpson','Evenness')), names_to = "Metric",values_to = "Alpha") %>% 
  mutate(Shrub = factor(Shrub, levels = shrublevs),
         Metric = factor(Metric, levels = c("Observed","Shannon","Simpson","Evenness")))

ggplot(gdat)+
  geom_boxplot(mapping=aes(x=Shrub, y=Alpha, fill=Shrub))+
  geom_jitter(mapping=aes(x=Shrub, y=Alpha), width = 0.2)+
  scale_fill_manual(values=shrubpal)+
  theme_bw()+
  facet_grid2(Metric~Site, scales='free', independent = T)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ggtitle('Alpha Diversity (ITS)')

ggsave('figs/vegan_shannon_panel_ITS.png', height=9, width=6)

######################################
ggplot(gdat)+
  geom_boxplot(mapping=aes(x=Shrub, y=Simpson, fill=Shrub))+
  geom_jitter(mapping=aes(x=Shrub, y=Simpson), width = 0.2)+
  theme_bw()+
  facet_wrap(~Site)+
  guides(fill='none')+
  ggtitle('Simpson Richness (vegan)')
ggsave('figs/vegan_simpson.png')
