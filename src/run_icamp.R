library(tidyverse)
library(iCAMP)
library(phyloseq)

#import 
ps <- readRDS('data/16Sdada2/phyloseq_dada2.tree.rds')

#wrangle that data
ptree <- phy_tree(ps)
asv_order <- ptree$tip.label
asv_table <- data.frame(otu_table(ps))[,asv_order]

print(identical(colnames(asv_table), ptree$tip.label))



### run the stegen et al framework

#generate phylogenetic - distance matrix
#pd_big <- pdist.big(tree = phy_tree(ps), 
#                    wd = 'data/icamp/', 
#                    nworker = 80,
#                    memory.G=200')
#pd <- cophenetic(ptree) #limited performance with high ASV 
#
#results <- qpen(comm = asv_table, 
#                pd = pd, 
#                tree = ptree,
#                nworker = 80,
#                memory.G=200,
#                rand = 1000) 
#
#saveRDS(results, "data/icamp/16S_icamp_qpen.rds")

###run with phylogenetic bins -- allows for taxa-specific effects
results <- icamp.big(comm = asv_table, 
                tree = ptree,
                pd.wd='data/icamp/',
                bin.size.limit=24,
                nworker = 80,
                memory.G=200,
                rand = 1000) 

saveRDS(results, "data/icamp/16S_icamp_bin.rds")

print(results$CbMPR)

