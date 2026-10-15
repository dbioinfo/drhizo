library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(colorspace)
library(vegan)
library(mirlyn)
setwd('~/WorkForaging/Academia/Nicole/drhizo')

#import data
ps <- t(readRDS('data/16Sdada2/phyloseq_dada2.rds'))
otu <- otu_table(ps, taxa_are_rows = F)
meta <- data.frame(sample_data(ps))
psdf <- phyloseq_to_df(t(ps))

#check out abundance plot to be sure
ibar <- bartax(ps, "Sample", taxrank="Phylum")
ibar

#rare curve
whole_rep <- rarefy_whole_rep(ps, rep = 100, mc.cores=10)
rarecurve(whole_rep, sample = "Sample")

#repeated subsampling rarefaction with mirl
mrare <- mirl(ps, libsize=2000, rep=100, replace=F, mc.cores=10)

#lets check the alpha diversity metrics first
alphadiv_df <- alphadivDF(mrare,diversity = "Shannon")
galpha <- alphawichVis(alphadiv_df, xvar = "SampleID", colorvar = "Shrub") + 
  theme(axis.text.x = element_text(angle=75, vjust = 0.7, hjust = .6))
galpha

acone <- alphacone(ps, rep=100, diversity = "shannon", replace=F, mc.cores=10)
gcone <- alphaconeVis(acone, cols = "Shrub")
gcone

#mirl beta diversity
betamatPCA_object <- betamatPCA(mrare, dsim = "bray", transformation = "hellinger")
betaPCA <- betamatPCA_object$x[,1:5]
betaPCA <- data.frame(PC1=betaPCA[,"PC1"], PC2=betaPCA[,"PC2"], PC3=betaPCA[,"PC3"], Sample=rownames(betaPCA))
betaPCA <- betaPCA %>% separate(Sample, into=c("rep","SampleID"), sep = '-')
betaPCA <- betaPCA %>% left_join(., meta, by="SampleID")
ggplot(betaPCA)+
  geom_point(aes(x=PC2, y=PC3, color=Shrub))+
  ggtitle("Bray-Curtis Beta Diversity PCA Rarefied")+
  theme_bw()

#bespoke beta diversity
#mirl only does PCA of beta diversity, so we'll need to improvise 
#their function just uses vegan under the hood anyways

btrans <- decostand(as.matrix(t(repotu_df(ps))), transformation="hellinger")

