library(tidyverse)
library(phyloseq)
library(ape)
library(Biostrings)
library(picante)
library(dada2)
setwd("~/WorkForaging/Academia/Nicole/drhizo/")

##write asv sequences to file
ps <- readRDS("data/16Sdada2/phyloseq_dada2.rds")
refs <- refseq(ps)
writeXStringSet(refs, "data/16Sqiime/asv-seqs.fasta", format = "fasta")

ps <- readRDS("data/ITSdada2/phyloseq_dada2.rds")
refs <- refseq(ps)
writeXStringSet(refs, "data/ITSqiime/asv-seqs.fasta", format = "fasta")

#run gen_newick

#import tree
ps <- readRDS("data/16Sdada2/phyloseq_dada2.rds")
ptree <- read.tree('data/16Sqiime/asv_ref_tree/tree.nwk')
asvs <- taxa_names(ps)
asvs <- intersect(asvs, ptree$tip.label)

ptree <- keep.tip(ptree, asvs)
ptree <- collapse.singles(ptree)

ps <- merge_phyloseq(ps, phy_tree(ptree))
saveRDS(ps,"data/16Sdada2/phyloseq_dada2.tree.rds")

#import tree
ps <- readRDS("data/ITSdada2/phyloseq_dada2.rds")
ptree <- read.tree('data/ITSqiime/asv_ref_tree/tree.nwk')
asvs <- taxa_names(ps)
asvs <- intersect(asvs, ptree$tip.label)

ptree <- keep.tip(ptree, asvs)
ptree <- collapse.singles(ptree)

ps <- merge_phyloseq(ps, phy_tree(ptree))
saveRDS(ps,"data/ITSdada2/phyloseq_dada2.tree.rds")


