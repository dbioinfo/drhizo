library(tidyverse)
library(dada2)
library(file2meco)
library(phyloseq)
library(microeco)
library(ComplexHeatmap)
library(circlize)
library(scales)  
setwd('~/WorkForaging/Academia/Nicole/drhizo')

ps <- readRDS('data/ITSdada2/phyloseq_dada2.rds')

otu <- otu_table(ps)

data <- phyloseq2meco(ps)
tf <- trans_func$new(data)
tf$cal_spe_func(prok_database = "FAPROTAX")
tf$cal_spe_func_perc(abundance_weighted = FALSE) 
saveRDS(tf, "data/faprotax/ITS_faprotax.rds")
head(tf$res_spe_func_perc)

#sort data into a complex heatmap
shrub_levels <- c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )
shrubpal <- c("#ED9F56","#7AC2CF","#A479A5","#5EBA49","#EFD273")

gdat <- tf$res_spe_func_perc
gdat <- gdat[which(colSums(gdat)>10)]
meta <- data.frame(sample_data(ps)) %>% mutate(Shrub = factor(Shrub, levels = shrub_levels))

mat <- t(as.matrix(gdat))
mat[is.na(mat)] <- 0  

normalize_rows_z <- function(mat) {
  t(apply(mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}
tmp <- normalize_rows_z(mat)

shrub_cols   <- setNames( shrubpal, shrub_levels)

ha <- HeatmapAnnotation(
  Shrub = meta$Shrub,
  col   = list(Shrub = shrub_cols),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(c(min(tmp), 0, max(tmp)), c("#1c005c","#f0f0f0", "#67000d"))
ht <- Heatmap(
  tmp[,meta %>% arrange(Shrub) %>% pull(SampleID)],
  name                 = "% FAPROTAX",
  column_split         = meta$Shrub,        # cluster within Shrub groups
  top_annotation       = ha,
  col                  = col_fun,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_column_names    = FALSE,
  row_title            = "",
  column_title         = "Samples",
  heatmap_legend_param = list(title = "Normalized Abundance")
)
ht
png("figs/ITSdada2/faprotax_heatmap.png", width = 2400, height = 1600, res = 300)
draw(ht, merge_legends = TRUE)
dev.off()

normalize_rows_z <- function(mat) {
  t(apply(mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}
