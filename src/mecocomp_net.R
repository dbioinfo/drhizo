library(phyloseq)
library(dplyr)
library(tidyverse)
library(vegan)
library(decontam)
library(file2meco)
library(SpiecEasi)
library(microeco)
library(meconetcomp)

setwd("~/WorkForaging/Academia/Nicole/")

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

#format for vegan
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

#filter data
conc <- colSums(asv)
contamdf <- isContaminant(ps, method="frequency",conc=conc, threshold=.1)
sum(contamdf$contaminant) 
ps_filt <- prune_taxa(!contamdf$contaminant, ps) #57 ASV removed at .1 thresh
ps_filt <- subset_samples(ps_filt, SampleID != 'SRR14874590' & SampleID != 'SRR14874591') #2 samples removed for low counts
ps_filt <- subset_taxa(ps_filt, Domain == "Bacteria") #filter for only Bacteria; 74 taxa removed 

#rarefy
otu <- as(otu_table(ps_filt), "matrix")
otu <- as.data.frame(otu)
rarecurve(otu, step = 100, cex = 0.75, las = 1, label = F) #???
conc #see lib sizes
rare <- rrarefy(otu, sample = 10000) #seeing some disagreement online by the maintainers of vegan that this may be unhelpful
tax <- as.matrix(tax)
rownames(tax) <- tax[,1]
ps_r <- phyloseq(otu_table(rare, taxa_are_rows = TRUE), subset_samples(mdat, SampleID %in% colnames(rare)), tax_table(tax) )

#construct networks separately
meco <- phyloseq2meco(ps_r)

shrubs <- unique(sample_data(ps_r)$Shrub)
netlist <- c()
for (shrub in shrubs){
  tmp <- clone(meco)
  tmp$sample_table <- tmp$sample_table %>% subset(Shrub == shrub)
  tmp$tidy_dataset()
  
  tnet <- trans_network$new(dataset=tmp,cor_method = 'spearman',filter_thres = 0.0005)
  tnet$cal_network(method='SpiecEasi')

  netlist[[shrub]] <- tnet
}

nsum <- cal_network_attr(netlist)
write.csv(nsum, 'data/16Snetworks/SpiecEasi_net_attr.csv', row.names = T)

netlist <- cal_module(netlist, undirected_method = 'cluster_fast_greedy')
netlist <- netlist %>% get_node_table() %>% get_edge_table()

#save the nets
saveRDS(netlist, file='data/16Snetworks/total_networks.rds')
netlist$Tarbush$save_network(filepath = 'data/16Snetworks/Tarbush_net.gexf')

#extract intersection
t1 <- edge_comp(netlist)
t2 <- trans_venn$new(t1)
t3 <- t2$trans_comm()
tintersect <- subset_network(netlist, venn=t3, name='Creosote&Mariola&Mesquite&Tarbush&Interspace')

t1 <- meco$merge_samples("Shrub")
t2 <- trans_venn$new(t1, ratio='seqratio')
p<-t2$plot_venn(petal_plot = TRUE)
ggsave('figs/ITS_SpiecEasi_venn_shrub.png', dpi=500, width=5000, height=4000, units='px')


