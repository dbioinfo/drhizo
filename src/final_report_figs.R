library(tidyverse)
library(ggh4x)
library(scales)
setwd('~/WorkForaging/Academia/Nicole/drhizo/')

shrubpal <- c("#ED9F56","#7AC2CF","#A479A5","#5EBA49","#EFD273")
shrublevs <- c("Creosote","Mariola","Mesquite","Tarbush","Interspace" ) 
meta <- read_csv("data/metadata.csv")
gdat <- meta %>% group_by(Site, Shrub) %>% summarize(n=n()) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

meta1 <- data.frame(sample_data(readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds"))) %>% rownames_to_column("Sample") %>% mutate(Amp="ITS")
meta2 <- data.frame(sample_data(readRDS("data/16Sdada2/phyloseq_dada2.OTU.rds"))) %>% rownames_to_column("SRA") %>% mutate(Amp="16S")
gdat <- rbind(meta1, meta2) %>% group_by(Site, Shrub, Amp) %>% summarize(n=n()) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

gg <- ggplot(gdat)+
  geom_point(aes(x=Shrub, y=Amp, size=n, color=n,shape=Amp))+
  geom_text(aes(x=Shrub, y=Amp, label=n))+
  scale_color_gradientn(colours = c('#DADADA',"#E7BC2E","#D85071"),values = c(0,.4,1 ) )+
  scale_size_continuous(range = c(5,8) )+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  guides(color='none', size='none',shape='none')+
  ylab("")+
  facet_grid2(Site ~ Shrub, scales = "free", switch = 'y',
              strip = strip_themed(background_x = list(element_rect(fill="white", colour = shrubpal[1], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[2], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[3], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[4], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[5], linewidth = 3)),
                                   background_y = list(element_rect(fill="white", colour = "#000000", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#000000", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#000000", linewidth = 2))))+
  ggtitle("Sample Sizes Post-QC")
gg  

ggsave("figs/sample_distribution.png", gg, width=7, height=3)



#####
## TAXONOMY Abundance
####
## 16S

ps <- readRDS("data/16Sdada2/phyloseq_dada2.OTU.rds")

ps_phy <- tax_glom(ps, taxrank = "Phylum")
ps_phy <- transform_sample_counts(ps_phy, function(x){x/sum(x)})
pdat <- psmelt(ps_phy) %>% 
  mutate(Phylum = case_when(Abundance < .05 ~"<5% Abundance",
                            .default = Phylum)) %>%  
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

randpal <- sample(RColorBrewer::brewer.pal(n = 10, name = "Paired"), 10, replace = F)

gg<-ggplot(pdat)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum, color=Kingdom),stat='identity',linewidth=.3)+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values =randpal)+
  scale_color_manual(values = c("darkorange","darkblue"))+
  facet_grid2(Shrub~Site, independent = "x", scales="free",
              strip = strip_themed(background_y = list(element_rect(fill="white", colour = shrubpal[1], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[2], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[3], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[4], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[5], linewidth = 3)),
                                   background_x = list(element_rect(fill="white", colour = "#00000080", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#00000080", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#00000080", linewidth = 2))))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(color = guide_legend(override.aes = list(fill=NA, linewidth=0.5)))
gg

ggsave("figs/final_abundance_plot.16S.png", gg, width=5, height=10)



## ITS

ps <- readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds")

ps_phy <- tax_glom(ps, taxrank = "Order")
ps_phy <- transform_sample_counts(ps_phy, function(x){x/sum(x)})
pdat <- psmelt(ps_phy) %>% 
  mutate(Order = case_when(Abundance < .10 ~"<10% Abundance",
                            .default = str_remove(Order, "o__") ),
         Phylum = case_when(Abundance < .10 ~"<10% Abundance",
                           .default = str_remove(Phylum, "p__"))  ) %>%  
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

randpal <- sample(RColorBrewer::brewer.pal(n = 12, name = "Paired"), 15, replace = T)

randpal <- sample(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#D85071", "#E31A1C",
             "#FDBF6F", "#EDAA56", "#FF7F00", "#CAB2D6", "#A479A5", "#6A3D9A", "#FFFF99",
             "#B15928"), 15)

minirandpal <- c("#1F78B4","#FF7F00", "#6A3D9A", "#E31A1C","#B15928")

gg<-ggplot(pdat)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Order, color=Phylum),stat='identity',linewidth=.3)+
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values = minirandpal)+
  scale_fill_manual(values =randpal)+
  facet_grid2(Shrub~Site, independent = "x", scales="free",
              strip = strip_themed(background_y = list(element_rect(fill="white", colour = shrubpal[1], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[2], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[3], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[4], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[5], linewidth = 3)),
                                   background_x = list(element_rect(fill="white", colour = "#00000080", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#00000080", linewidth = 2), 
                                                       element_rect(fill="white", colour = "#00000080", linewidth = 2))))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(color = guide_legend(override.aes = list(fill=NA, linewidth=0.5)))
gg

ggsave("figs/final_abundance_plot.ITS.png", gg, width=5, height=10)




## specific taxa plots
#16S 
gg<- ggplot(pdat )+
  geom_boxplot(aes(x=Phylum, y=Abundance, fill=Phylum), outliers = F)+
  geom_jitter(aes(x=Phylum, y=Abundance), height=0, width=0.3, alpha=0.6, size=1)+
  scale_fill_manual(values =randpal)+
  scale_y_continuous(expand = c(0, 0))+
  facet_grid2(~Shrub, independent = "x", scales="free",
              strip = strip_themed(background_x = list(element_rect(fill="white", colour = shrubpal[1], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[2], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[3], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[4], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[5], linewidth = 3))))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")
gg
ggsave("figs/phylum_abundance_boxplot.16S.png", height=5, width=10)


# ITS
pdat <- pdat %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
gg<- ggplot(pdat )+
  geom_boxplot(aes(x=Order, y=Abundance, fill=Order), outliers = F)+
  geom_jitter(aes(x=Order, y=Abundance), height=0, width=0.3, alpha=0.6, size=1)+
  scale_fill_manual(values =randpal)+
  scale_y_continuous(expand = c(0, 0))+
  facet_grid2(~Shrub, independent = "x", scales="free",
              strip = strip_themed(background_x = list(element_rect(fill="white", colour = shrubpal[1], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[2], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[3], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[4], linewidth = 3),
                                                       element_rect(fill="white",colour = shrubpal[5], linewidth = 3))))+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")
gg
ggsave("figs/order_abundance_boxplot.ITS.png", height=5, width=10)




### faprotax 
##16S
ps <- readRDS("data/16Sdada2/phyloseq_dada2.OTU.rds")
tf  <- readRDS('data/faprotax/16S_faprotax.rds')
gdat <- tf$res_spe_func_perc
gdat <- gdat[which(colSums(gdat)>10)]
meta <- data.frame(sample_data(ps)) %>% mutate(Shrub = factor(Shrub, levels = shrub_levels)) %>% arrange(Shrub)

mat <- t(as.matrix(gdat))[,meta %>% rownames()]
mat[is.na(mat)] <- 0  

normalize_rows_z <- function(mat) {
  t(apply(mat, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
}
tmp <- normalize_rows_z(mat)

shrub_cols   <- setNames( shrubpal, shrub_levels)

ha <- HeatmapAnnotation(
  Shrub = meta$Shrub,
  col   = list(Shrub = shrub_cols),
  annotation_name_side = "left",
  annotation_legend_param = list(Shrub = list(at = shrub_levels))
)

col_fun <- colorRamp2(c(-3, 0, 3), c("#1c005c","#f0f0f0", "#67000d"))
ht <- Heatmap(
  tmp,
  name                 = "% FAPROTAX",
  column_split         = meta$Shrub,        # cluster within Shrub groups
  top_annotation       = ha,
  col                  = col_fun,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_column_names    = FALSE,
  row_title            = "",
  column_title         = "Samples", cluster_columns = F,
  heatmap_legend_param = list(title = "Normalized Abundance")
)
ht
png("figs/16Sdada2/faprotax_heatmap.png", width = 2400, height = 1600, res = 300)
draw(ht, merge_legends = TRUE)
dev.off()



###ITS2
ps <- readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds")
tf  <- readRDS('data/faprotax/ITS_faprotax.rds')
gdat <- tf$res_spe_func_perc
gdat <- gdat[which(colSums(gdat)>10)]
meta <- data.frame(sample_data(ps)) %>% mutate(Shrub = factor(Shrub, levels = shrub_levels)) %>% arrange(Shrub)

mat <- t(as.matrix(gdat))[,meta %>% arrange(Shrub) %>% pull(Sample)]
mat[is.na(mat)] <- 0  

tmp <- normalize_rows_z(mat)

shrub_cols   <- setNames( shrubpal, shrub_levels)

ha <- HeatmapAnnotation(
  Shrub = meta$Shrub,
  col   = list(Shrub = shrub_cols),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(c(-3, 0, 3), c("#1c005c","#f0f0f0", "#67000d"))
ht <- Heatmap(
  tmp,
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



### beta diversity 
ps <- readRDS("data/16Sdada2/phyloseq_dada2.tree.rds")
meta <- data.frame(sample_data(ps)) 
ord <- ordinate(ps, method = "PCoA", distance = "unifrac")
plot_ordination(ps, ord, color="Shrub") + geom_point(size=3) +scale_color_manual(values = shrubpal)

gdat <- ord$vectors[,1:3] 
colnames(gdat) <- c("PCoA1","PCoA2","PCoA3")
gdat <- data.frame(gdat) %>% rownames_to_column("SampleID") %>% left_join(., meta, by="SampleID") %>% 
  mutate(Shrub=factor(Shrub, levels=shrublevs))
centroids <- gdat %>% group_by(Shrub) %>% summarise(across(all_of(c('PCoA1', 'PCoA2', 'PCoA3')), function(x){mean(x)}) )
gg <- ggplot(gdat)+
  geom_point(aes(x=PCoA1, y=PCoA2, color = Shrub, shape = Site), size=4)+
  stat_ellipse(aes(x=PCoA1, y=PCoA2, color = Shrub))+
  geom_point(data=centroids, aes(x=PCoA1, y=PCoA2, color = Shrub), size=6, shape=4, stroke=2)+
  scale_color_manual(values = shrubpal)+
  theme_bw()+
  ggtitle("Unifrac PCoA (16S Unweighted)")
gg
ggsave("figs/16S_PCoA_Unifrac.png",gg)



ps <- readRDS("data/ITSdada2/phyloseq_dada2.tree.rds")
meta <- data.frame(sample_data(ps)) 
ord <- ordinate(ps, method = "PCoA", distance = "unifrac")

gdat <- ord$vectors[,1:3] 
colnames(gdat) <- c("PCoA1","PCoA2","PCoA3")
gdat <- data.frame(gdat) %>% rownames_to_column("SampleID") %>% left_join(., meta, by="SampleID") %>% 
  mutate(Shrub=factor(Shrub, levels=shrublevs))
centroids <- gdat %>% group_by(Shrub) %>% summarise(across(all_of(c('PCoA1', 'PCoA2', 'PCoA3')), function(x){mean(x)}) )
gg <- ggplot(gdat)+
  geom_point(aes(x=PCoA1, y=PCoA2, color = Shrub, shape = Site), size=4)+
  stat_ellipse(aes(x=PCoA1, y=PCoA2, color = Shrub))+
  geom_point(data=centroids, aes(x=PCoA1, y=PCoA2, color = Shrub), size=6, shape=4, stroke=2)+
  scale_color_manual(values = shrubpal)+
  theme_bw()+
  ggtitle("Unifrac PCoA (ITS Unweighted)")
gg
ggsave("figs/ITS_PCoA_Unifrac.png",gg)





##Spiec Easi figs
dat16 <- read_csv("data/networks/SpiecEasi_union_attr.csv") %>% 
  rename("Module"="stat_inf_class") %>% mutate(Amp="16S") %>% 
  select(Label, Phylum=att2, Cum_Abundance=att3, eigencentrality, clustering, Module, Amp)
datit <- read_csv("data/ITSnetworks/SpiecEasi_union_attr.csv") %>% 
  rename("Module"="att4") %>% mutate(Amp="ITS") %>% 
  select(Label, Phylum=att2, Cum_Abundance=att3, eigencentrality, clustering, Module, Amp)

gdat <- rbind(dat16, datit)

tmp <- gdat %>% group_by(Amp, Module) %>% summarise(n=n()) %>% 
  filter( ((Amp=="16S")&(n>=5))|((Amp=="ITS")&(n>=3)) ) %>% 
  arrange(-n)
xord <- tmp$Module

g1 <- ggplot(tmp )+
  geom_bar(aes(x=reorder(Module, -n), y=n, fill=Amp), stat='identity', )+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(~Amp, nrow=1, scales='free_x')+
  theme_bw()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill= 'none' )+
  ggtitle("Shrub Network Module Sizes")+
  xlab("Module")+
  ylab("N")
g1

tmp <- gdat %>% filter( ((Amp=="16S")&(eigencentrality>=.9))|((Amp=="ITS")&(eigencentrality>=.8)) )

g2 <- ggplot(tmp %>% filter(Amp=="16S")) +
  geom_point(aes(x=eigencentrality, y=clustering,color=Phylum,size=Cum_Abundance))+
  scale_color_manual(values = sample(pal_brewer( palette = "Spectral", direction = 1)(8), 8, replace = F))+
  scale_size_continuous(range = c(2,5), name="Cumulative Relative Abundance")+
  theme_bw()+
  ggtitle("16S Shrub Network Centrality")+
  xlab("Eigencentrality")+ylab("Clustering Coefficient")
g2

g3 <- ggplot(tmp %>% filter(Amp=="ITS")) +
  geom_point(aes(x=eigencentrality, y=clustering,color=Phylum,size=Cum_Abundance))+
  scale_color_manual(values = sample(pal_brewer( palette = "Spectral")(8), 3, replace = F))+
  scale_size_continuous(range = c(2,5), name="Cumulative Relative Abundance")+
  theme_bw()+
  ggtitle("ITS Shrub Network Centrality")+
  xlab("Eigencentrality")+ylab("Clustering Coefficient")
g3

library(patchwork)
png("figs/netfigs/shrub_bottom_panel_stats.png", width = 1000, height=500, units = "px")
(g1 +g2 +g3)
dev.off()



ttmp <- readRDS("data/16SIndSpec/res.shrub.interspace.16S.rds")
