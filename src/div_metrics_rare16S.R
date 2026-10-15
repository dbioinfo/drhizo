library(tidyverse)
library(dada2)
library(phyloseq)
library(vegan)
setwd('~/WorkForaging/Academia/Nicole/drhizo')

#import some stuff
ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')
meta <- sample_data(ps)
otu <- data.frame(otu_table(ps))

#rarefaction curve
rlev <- 2500
gg <- rarecurve(x = otu, step=10, tidy = T)
ggplot(gg %>% group_by(Site)) +
  geom_line(aes(x=Sample, y=Species, group=Site))+
  geom_vline(mapping = aes(xintercept=rlev))

smry <- gg %>% group_by(Site) %>% summarize(msamp = max(Sample))
minlib <- min(smry$msamp) #smallest library size
niter <- round(minlib/rlev,0) #number of without-replacement iterations possible before removing a sample
summary(smry$msamp)

#rarefaction alpha
r1 <- rarefy(otu, rlev)
for (i in 1:niter){
  
  
}


#rarefaction beta

#rarefaction ancom-bc