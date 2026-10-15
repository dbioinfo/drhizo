library(tidyverse)
library(dada2)
library(file2meco)
library(phyloseq)
library(microeco)
library(ComplexHeatmap)
library(circlize)
library(scales)  
setwd('~/WorkForaging/Academia/Nicole/drhizo')

ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')

otu <- otu_table(ps)

data <- phyloseq2meco(ps)
tf <- trans_func$new(data)
tf$cal_spe_func(prok_database = "FAPROTAX")
tf$cal_spe_func_perc(abundance_weighted = FALSE) 

head(tf$res_spe_func_perc)

#sort data into a complex heatmap
gdat <- tf$res_spe_func_perc
gdat <- gdat[which(colSums(gdat)>10)]
meta <- sample_data(ps)

mat <- t(as.matrix(gdat))
mat[is.na(mat)] <- 0  

shrub_levels <- unique(meta$Shrub)
shrub_cols   <- setNames(hue_pal()(length(shrub_levels)), shrub_levels)

ha <- HeatmapAnnotation(
  Shrub = meta$Shrub,
  col   = list(Shrub = shrub_cols),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(c(0, max(mat)), c("#f0f0f0", "#67000d"))
ht <- Heatmap(
  mat,
  name                 = "% FAPROTAX",
  column_split         = meta$Shrub,        # cluster within Shrub groups
  top_annotation       = ha,
  col                  = col_fun,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_column_names    = FALSE,
  row_title            = "",
  column_title         = "Samples",
  heatmap_legend_param = list(title = "Relative Abundance")
)
png("figs/faprotax_heatmap.png", width = 2400, height = 1600, res = 300)
draw(ht, merge_legends = TRUE)
dev.off()

