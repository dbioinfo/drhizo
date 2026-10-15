library(tidyverse)
library(dada2)
library(file2meco)
library(phyloseq)
library(microeco)
library(scales)  

setwd('~/WorkForaging/Academia/Nicole/drhizo/')

ps <- readRDS('data/16Sdada2/phyloseq_dada2.tree.rds')
otu <- data.frame(otu_table(ps))
meta <- data.frame(sample_data(ps))

data <- phyloseq2meco(ps)
tf <- trans_func$new(data)
tf$cal_spe_func(prok_database = "FAPROTAX")
tf$cal_spe_func_perc(abundance_weighted = FALSE) 


head(tf$res_spe_func_perc)
saveRDS(tf, 'data/faprotax/16S_faprotax.rds')
