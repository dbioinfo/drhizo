library(tidyverse)
library(colorspace)
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
meta <- read_csv('data/SraRunTable.csv') %>% rename(SRAid = "Run")

#Phylum level summary
tdata <- data %>% 
  pivot_longer(-OTU, names_to = "SRAid", values_to = "count") %>% 
  filter(count>0) %>% 
  left_join(tax %>% select(OTU, Phylum), by="OTU") %>% 
  group_by(SRAid, Phylum) %>% 
  summarize(count=sum(count)) %>%
  group_by(SRAid) %>%  
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>% 
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) 

#select > 5% Phyla
rteps <- tdata %>% mutate(Phylum=case_when(percent<1~"<1% Abundance",
                                           .default = Phylum)) %>% 
  ungroup() %>% 
  group_by(SRAid, Phylum) %>% 
  summarize(count=sum(count)) %>%
  ungroup() %>% 
  group_by(SRAid) %>% 
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>%   
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) 

#add metadata
rteps <- rteps %>% left_join(meta %>% select(Site, Shrub, SRAid, Sample_ID), by='SRAid') 

#set up graph and save
palette <- qualitative_hcl(18, palette = "Dynamic") # or "Set3", "Classic Dark"
plevs <- unique(rteps$Phylum)
palette <- setNames(palette, nm = plevs)

ggplot(rteps) +
  geom_bar(aes(x=Sample_ID, y=percent, fill=Phylum), stat='identity', na.rm = T) +
  facet_wrap(~Site+Shrub, nrow=1, scales = "free_x") +
  scale_fill_manual(values = palette)+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("Sample")+
  ylab("Relative Abundance")

ggsave("figs/Phylum_Rel_Abundance.png", dpi=300, width=5000, height=2000, units='px')  
