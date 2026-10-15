library(tidyverse)
library(iCAMP)
library(phyloseq)

setwd('~/WorkForaging/Academia/Nicole/drhizo/')

res <- readRDS('data/icamp/ITS/ITS_icamp_bin.rds')
ps <- readRDS('data/ITSdada2/phyloseq_dada2.tree.rds')
meta <- data.frame(sample_data(ps))


#global patterns
global_proc <- res$CbMPDiCBraya %>% 
  pivot_longer(cols = c(Heterogeneous.Selection, Homogeneous.Selection, 
                        Dispersal.Limitation, Homogenizing.Dispersal, 
                        Drift.and.Others),
               names_to = "Process",
               values_to = "Relative_Importance") %>% 
  rename(SampleID=sample1) %>%
  left_join(., meta, by="SampleID") %>% 
  group_by(Shrub, Process) %>% 
  summarize(Relative_Importance=mean(Relative_Importance)) %>% 
  mutate(Shrub=factor(Shrub, levels=c("Creosote","Mariola","Tarbush","Mesquite","Interspace")))

ggplot(global_proc) +
  geom_bar(mapping=aes(x=Shrub, y=Relative_Importance, fill=Process), stat='identity')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = c("#B35C1E","#178764","#5C1EB3","#A61EB3","#93AD1F"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  xlab("")+ ylab("Relative Importance")+
  ggtitle("Global Patterns of Community Assembly (16S)")
ggsave('figs/icamp/16S/global_CA_16S.png', width=7, height=5)

global_proc %>% group_by(Process) %>% summarize(meanperc = mean(Relative_Importance))



#across-bin analysis
weights <- res$detail$bin.weight
sig_bnti <- res$detail$SigbMPDi
sig_rc <- res$detail$SigBCa

tmp <- res$CbMPDiCBraya %>% 
  pivot_longer(cols = c(Heterogeneous.Selection, Homogeneous.Selection, 
                        Dispersal.Limitation, Homogenizing.Dispersal, 
                        Drift.and.Others),
               names_to = "Process",
               values_to = "Relative_Importance") %>% 
  group_by(sample1, sample2) %>% 
  slice_max(order_by = Relative_Importance, n = 1) 


#calculate abundance-weighted metrics across phylo bins
wBNTI <- rowSums(sig_bnti[, -c(1,2)] * weights[, -c(1,2)], na.rm = TRUE)
wRC <- rowSums(sig_rc[, -c(1,2)] * weights[, -c(1,2)], na.rm = TRUE)
plot_df <- data.frame(sample1=sig_bnti$name1, sample2=sig_bnti$name2,BNTI = wBNTI,RC = wRC) %>% 
  left_join(., tmp, by=c('sample1','sample2'))
ggplot(plot_df %>% arrange(Process))+
  geom_point(mapping=aes(x=BNTI, y=RC, color=Process), alpha=0.8)+
  theme_bw()+
  scale_x_continuous(limits = c(-1,1))+
  scale_y_continuous(limits = c(-1,1))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  scale_color_manual(values = c("#B35C1E","#178764","#5C1EB3"))+
  ggtitle("Process Dominance (16S Sample-Sample)")
ggsave('figs/icamp/16S/process_dominance_16S.png', width = 7, height=5)


#investigate pairs with non-drift 
plot_df <- plot_df %>% 
  left_join(., meta %>% select(sample1=SampleID, Shrub.S1=Shrub), by='sample1') %>% 
  left_join(., meta %>% select(sample2=SampleID, Shrub.S2=Shrub), by='sample2') 
  
##homogenous selection
HOS <- plot_df %>% filter(Process=='Homogeneous.Selection')
table(HOS$Shrub.S1, HOS$Shrub.S2)

##dispersal limitation
DL <- plot_df %>% filter(Process=='Dispersal.Limitation')
table(DL$Shrub.S1, DL$Shrub.S2)

#which ASVs drive those patterns
asv_to_bin <- res$detail$taxabin$sp.bin %>% rownames_to_column("ASV")
tax <- data.frame(tax_table(ps) )%>% rownames_to_column("ASV")
phylo_bin <- left_join(asv_to_bin, tax, by="ASV")  %>% rename(bin=bin.id.new) %>% mutate(bin=paste0('bin',bin))

phy_bnti <- sig_bnti %>% 
  pivot_longer(!c(name1, name2), names_to = 'bin', values_to = 'CbMPDi') %>% 
  mutate(bin=sub('CbMPDi.','',bin))
phy_rc <- sig_rc  %>% 
  pivot_longer(!c(name1, name2), names_to = 'bin', values_to = 'Cbraya') %>% 
  mutate(bin=sub('Cbraya.','',bin))


#filter to phyla with mean >0.005 rel abundance
ps_tmp <- transform_sample_counts(ps, function(x) {x / sum(x)})
ps_tmp <- tax_glom(ps_tmp, taxrank = 'Order')
ps_tmp <- filter_taxa(ps_tmp, function(x) mean(x) > 0.005, TRUE)
pnames <- as.character(tax_table(ps_tmp)[,"Order"])

#one samp comparison from HOS
gdat <- left_join(phy_bnti, phy_rc, by=c('name1','name2','bin'))
gdat <- left_join(gdat, phylo_bin, by='bin')
tdat <- gdat %>% filter(name1=="SRR14874576",name2=="SRR14874551",
                        Order %in% pnames) %>% 
  group_by(name1, name2, bin, Order) %>% 
  summarize(CbMPDi=mean(CbMPDi), Cbraya=mean(Cbraya)) %>% 
  mutate(Order=str_remove(Order, pattern="o__"))
#filter(name1=="M49",name2=="M34",

#set up background area to help comprehend graph
rectdat <- bind_rows(list(xmin=-.95, ymin=-.95,xmax=.95,ymax=.95, DominantProcess='Drift'),
                       list(xmin=-Inf, ymin=-Inf,xmax=-.95,ymax=Inf, DominantProcess='Selection'),
                       list(xmin=.95, ymin=-Inf,xmax=Inf,ymax=Inf, DominantProcess='Selection'),
                       list(xmin=-.95, ymin=.95,xmax=.95,ymax=Inf, DominantProcess='Dispersal'),
                       list(xmin=-.95, ymin=-Inf,xmax=.95,ymax=-.95, DominantProcess='Dispersal') )
randpal <- sample(c("#8A59EA","#DC59EA","#F05BB1",
                    "#E45779","#E46957","#CC513C",
                    "#F3B166","#D9F366","#A2F366",
                    "#54DFC8","#43B26A","#439AB2","#5971EA",
                    "#5002A4", "#7F0089", "#8C015F", "#841136" ,"#802711", 
                    "#69271C", "#7D5202", "#647400", "#417500", "#006C5F",
                    "#17582F", "#074D5C", "#00319F"))
#graph single sample comparison
ggplot(tdat )+
  geom_vline(xintercept = c(-.95,.95), color='black', linewidth=1.5)+
  geom_segment(x=-.95, y=.95, xend=.95, yend=.95, color='black', linewidth=1.5)+
  geom_segment(x=-.95, y=-.95, xend=.95, yend=-.95, color='black', linewidth=1.5)+
  geom_rect(data=rectdat, mapping=aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=DominantProcess),
            alpha=0.25)+
  scale_fill_manual(values = rev(c('lightgreen',"grey",'lightblue')), name="Dominant Process" )+
  geom_jitter(mapping=aes(x=CbMPDi, y=Cbraya, color=Phylum), size=3, alpha=0.7, height=0.05, width=0.05)+
  scale_color_manual(values= randpal)+
  theme_bw()+
  theme(axis.text = element_text(color='black'),
        axis.title = element_text(face='bold',color='black'))+
  guides(color=guide_legend(ncol=2))+
  ggtitle("Phylum-level Effects", subtitle = "ITS |  M49 (Interspace) vs M34 (Creosote)")

ggsave('figs/icamp/16S/M49.M34.Phylum.png', width=7, height=5)


#for homogeneous selection, we look for which phylum dominates the bottom left quadrant. 
#this is because in homogeneous selection species are more closely related between samples than expected by chance (negative CbMPDi)
#and also, their relative abundance is very similar (negative Cbraya and above 0.005)
#here, the two samples show that Phylum Actinomycetota dominates that quadrant, and is therefore likely being selected for in both samples
#they're both Mariola samples, so there may be a relationship between that shrub and Actinomycetota


#we can do a similar analysis between two microhabitats (Mariola vs others)
tdat <- gdat %>% 
       left_join(., meta %>% select(name1=SampleID, Shrub.S1=Shrub), by='name1') %>% 
       left_join(., meta %>% select(name2=SampleID, Shrub.S2=Shrub), by='name2')  %>% 
       filter(Order %in% pnames)%>% 
       mutate(Order=str_remove(Order, pattern="o__"))
t1 <- tdat %>% filter(Shrub.S1=="Interspace",Shrub.S2!="Interspace")
t2 <- tdat %>% filter(Shrub.S2=="Interspace",Shrub.S1!="Interspace")
tdat <- rbind(t1, t2) %>% 
  group_by(Shrub.S1,Shrub.S2, bin, Order) %>% 
  summarize(CbMPDi=mean(CbMPDi),Cbraya=mean(Cbraya)) %>% 
  mutate(Contrast = case_when(Shrub.S1=="Interspace"~Shrub.S2,
                              Shrub.S2=="Interspace"~Shrub.S1,
                              .default=NA)) %>% 
  mutate(Contrast = factor(Contrast, levels=c("Creosote","Mariola","Tarbush","Mesquite")))

randpal <- sample(c("#8A59EA","#DC59EA","#F05BB1",
                    "#E45779","#E46957","#CC513C",
                    "#F3B166","#D9F366","#A2F366",
                    "#54DFC8","#43B26A","#439AB2","#5971EA",
                    "#5002A4", "#7F0089", "#8C015F", "#841136" ,"#802711", 
                    "#69271C", "#7D5202", "#647400", "#417500", "#006C5F",
                    "#17582F", "#074D5C", "#00319F"))

ggplot(tdat)+  
  geom_vline(xintercept = c(-.95,.95), color='black', linewidth=1.5)+
  geom_segment(x=-.95, y=.95, xend=.95, yend=.95, color='black', linewidth=1.5)+
  geom_segment(x=-.95, y=-.95, xend=.95, yend=-.95, color='black', linewidth=1.5)+
  geom_rect(data=rectdat, mapping=aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=DominantProcess),
            alpha=0.25)+
  scale_fill_manual(values = rev(c('lightgreen',"grey",'lightblue')), name="Dominant Process" )+
  geom_jitter(mapping=aes(x=CbMPDi, y=Cbraya, color=Order,shape=Contrast), size=3, alpha=0.6, height=0.05, width=0.05)+
  theme_bw()+
  scale_color_manual(values= randpal)+
  guides(color=guide_legend(order=2, ncol=2))+
  theme(legend.box = "vertical")+
  ggtitle("Order-level Community Assembly", subtitle = "ITS | Interspace vs Shrubs")

ggsave('figs/icamp/16S/Interspace_vs_all_CA_Order.png', width = 9, height=7)


#specific functions of bins in HeS
#tdat is interspace vs others
funcdat <- readRDS('data/faprotax/16S_faprotax.rds')
asv_func <- funcdat$res_spe_func %>% rownames_to_column("ASV")
asv_func <- asv_func %>% 
  left_join(., asv_to_bin %>% select(ASV, bin=bin.id.new), by="ASV") %>% 
  pivot_longer(!c(ASV, bin), names_to = 'Function', values_to = 'n_asv') %>% 
  group_by(bin, Function) %>% 
  summarize(n_asv = sum(n_asv)) %>% 
  filter(n_asv>0) %>% 
  mutate(bin=paste0('bin',bin)) %>% 
  left_join(., gdat %>% select(bin, Phylum), by='bin')

asv_func <- asv_func  %>% unique() %>% mutate(Phylum  = str_remove(Phylum, "o__")) 

hes_bins <- tdat %>% filter(CbMPDi> 0.5) %>% pull(bin) %>% unique()
hos_bins <- tdat %>% filter(CbMPDi< -0.5) %>% pull(bin) %>% unique()


ggplot(asv_func %>% filter(bin %in% hes_bins))+
  geom_tile(mapping=aes(x=bin, y=Function, fill=n_asv))+
  scale_fill_viridis_c(name = "N ASVs", option = "C", limits = c(0,30), oob=scales::squish)+
  scale_y_discrete(limits=rev)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  xlab("") + ylab("")+
  ggtitle("Functional Annotation Heatmap",subtitle = "ITS | Interspace Heterogeneous Selected")

ggsave("figs/icamp/ITS/func.interspace.HeS.png", width=7, height=7)

otu <- data.frame(t(otu_table(ps %>% transform_sample_counts(fun = function(x){x / sum(x)}))) ) %>% rownames_to_column("ASV")
heatdat <- otu %>% 
  left_join(., asv_to_bin %>% mutate(bin=paste0('bin',bin.id.new) ) %>% select(ASV, bin) , by="ASV") %>% 
  filter(bin %in% hes_bins) %>% 
  group_by(bin) %>% summarize(across(starts_with("M"), mean)) %>% 
  pivot_longer(!bin, names_to = "SampleID", values_to = "Rel_Abundance") %>% 
  left_join(., meta, by="SampleID") %>% 
  ungroup() %>% group_by(Shrub, bin) %>% 
  summarize(Mean_Rel_Abundance=mean(Rel_Abundance)) %>% 
  left_join(., asv_func %>% select(bin, Function) %>% unique(), by="bin") %>% 
  mutate(Shrub = factor(Shrub, levels=c("Creosote","Mariola","Tarbush","Mesquite","Interspace")))

ggplot(heatdat)+
  geom_tile(mapping=aes(x=bin, y=Shrub, fill=log10(Mean_Rel_Abundance)))+
  scale_fill_viridis_c(name = "log10(Mean RA)", option = "C", 
                       limits = c(-4,-3), oob=scales::squish)+
  scale_y_discrete(expand=c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text())+
  xlab("")+ylab("")+
  ggtitle("Phylogenetic Bin Abundance",subtitle = "ITS | Interspace Heterogeneous Selected")
ggsave("figs/icamp/ITS/phybin.abund.interspace.HeS.png", width=7, height=7)


#are bins monophyletic? no. barely ever. 
library(ape)
ptree <- phy_tree(ps)
for (ibin in unique(asv_to_bin$bin)){

  iasvs <- asv_to_bin %>% filter(bin==ibin) %>% pull(ASV) %>% unique()
  print(paste0(ibin, " monophyly: ", is.monophyletic(ptree, iasvs)))
  
  
}


