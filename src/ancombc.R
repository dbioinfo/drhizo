library(tidyverse)
library(colorspace)
library(phyloseq)
library(ANCOMBC)
setwd("WorkForaging/Academia/Nicole/")

#import data
tax <- read_tsv('data/ASV/tax_dir/taxonomy.tsv')  %>% rename(OTU="Feature ID")  %>%
  separate(
    col = Taxon,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    fill = "right",   # if fewer than 7 levels, fill with NA
    remove = FALSE    
  ) %>%
  mutate(across(Domain:Species, ~ sub("^[a-z]__","",.)))
data <- read_tsv('data/ASV/feature_table.tsv', skip=1) %>% rename(OTU="#OTU ID")
meta <- read_csv('data/SraRunTable.csv') %>% rename(SampleID = "Run")

#format for ancombc
otu_id <- data$OTU
asv <- otu_table(as.matrix(data[, 2:ncol(data)]), taxa_are_rows = T)
rownames(asv) <- data$OTU
mdat <- sample_data(meta %>% 
                      select(SampleID, Site, Shrub) %>% 
                      mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush"))))
rownames(mdat) <- mdat$SampleID
ttable <- tax_table(as.matrix(tax %>% select(-c("OTU","Taxon","Confidence"))))
rownames(ttable) <- tax$OTU
ps <- phyloseq(asv, mdat, ttable)

#run it
out <- ancombc2(data = ps,
               fix_formula='Site+Shrub',
               rand_formula = NULL,
               tax_level="Genus",
               p_adj_method = 'holm',
               prv_cut = 0.1,
               n_cl = 10,
               group='Shrub',
               pairwise = T,
               struc_zero = T,
               neg_lb = F)

#graph results
out$res_pair %>% select(starts_with("diff")) %>% summarize_all(sum)
gdat <- out$res_pair
sigtaxs<-gdat %>% filter(diff_ShrubCreosote|diff_ShrubMariola|diff_ShrubMesquite|diff_ShrubTarbush) %>% pull(taxon)
ggdat <- gdat %>% 
  filter(taxon %in% sigtaxs) %>% 
  group_by(taxon) %>% 
  select(starts_with("lfc")) %>% 
  pivot_longer(starts_with("lfc"),names_to='comparison', values_to = 'LFC') %>% 
  mutate(comparison = sub("^lfc_", "", comparison)) %>% 
  mutate(comparison= sub("*Shrub*","", comparison))%>% 
  mutate(comparison= sub("*Shrub*","", comparison))

glevs <- c("Creosote","Mariola", "Mesquite","Tarbush", "Mariola_Creosote", 
           "Mesquite_Creosote", "Tarbush_Creosote",  "Mesquite_Mariola","Tarbush_Mariola","Tarbush_Mesquite" )

ggdat <- ggdat %>% mutate(comparison=factor(comparison, levels=glevs))

gg <- ggplot(ggdat)+
  geom_tile(mapping=aes(x=comparison, y=taxon, fill=LFC))+
  scale_fill_viridis_c(option='plasma')+
  theme_bw() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(), 
    panel.border = element_blank())+
  xlab("")+
  ylab("Taxon")
gg
ggsave('figs/Genus_lvl_heatmap.png', height=4000, width=4000, dpi=300, units='px')
