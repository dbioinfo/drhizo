library(tidyverse)
library(phyloseq)
library(indicspecies)
library(permute)
library(ggrastr)
library(Hmisc)
library(ggh4x)


setwd('~/WorkForaging/Academia/Nicole/drhizo/')

shrubpal <- c("#ED9F56","#7AC2CF","#A479A5","#5EBA49","#EFD273")
shrublevs <- c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )

###
# 16S 
###

#import
ps <- readRDS("data/16Sdada2/phyloseq_dada2.OTU.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- t(data.frame(otu_table(ps)))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV") %>% mutate(across(Kingdom:Species, function(x){str_remove(x, '[a-z]__')}))

#ID unique taxa, remove them
uasv <- (asvtab > 0)
uasv <- names(which(colSums(uasv)>1))

asvtab <- data.frame(asvtab) %>% select(all_of(uasv))

#run program
isares.shrub <- multipatt(asvtab, meta$Shrub,
                          control=how(nperm = 1000), max.order = 2)
saveRDS(isares.shrub,"data/16SIndSpec/res.shrub.16S.rds")
#isares.shrub <- readRDS("data/16SIndSpec/res.shrub.rds")

#interpret
spec <- data.frame(isares.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity")

sens <- data.frame(isares.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")

scores <- data.frame(isares.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")

gdat <- left_join(spec, sens, by=c("ASV","Shrub") ) %>% left_join(., scores, by=c("ASV","Shrub"))

tdat <- gdat %>% filter(Shrub %in% shrublevs ) %>% 
  mutate(Shrub = factor(Shrub, levels = shrublevs))
g1a <- ggplot(tdat)+
  geom_function(fun=function(x){(0.75**2)/x}, color='red')+
  geom_point_rast(aes(x=specificity, y=sensitivity, fill=Shrub),alpha=0.8, shape=21, stroke=0.25,size=3)+
  scale_y_continuous(expand=c(0,0), limits = c(0,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=shrubpal)+
  theme_bw()+
  xlab("Specificity")+
  ylab("Sensitivity")+
  ggtitle("Indicator Species Analysis (Microhabitat Single Level 16S)")
g1a
ggsave("figs/isa_shrub_simple.16S.png", g1a)

#now with phi statistic
phires.shrub <- multipatt(asvtab, meta$Shrub, func = 'r.g',
                          control=how(nperm = 1000), max.order = 2)
saveRDS(phires.shrub,"data/16SIndSpec/res.phi.shrub.16S.rds")

pscores <- as.data.frame(phires.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

g2a <- ggplot(pscores %>% filter(Shrub %in% c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )))+
  geom_jitter(aes(x = Shrub, y=Phi, fill=abs(Phi), alpha=abs(Phi) ), 
              height=0, shape=21,size=2)+
  scale_fill_gradientn(colors = c("black","darkred","darkorange","yellow"), name = "| Phi |")+
  scale_alpha_continuous(range =  c(0,0.9))+
  guides(alpha='none')+
  ggtitle("Phi Correlation (16S OTUs across Shrubs)")+
  theme_bw()
g2a

ggsave("figs/phi_shrub_simple_16S.png", g2a)



### what's with those horizontal lines? -- they were asvs unique to a sample
linedat <- pscores %>% filter(Shrub %in% c("Creosote","Mariola","Mesquite","Tarbush","Interspace" ),Phi>.25 ) %>% 
  group_by(Shrub, Phi) %>% summarise(n=n()) %>% 
  filter(n>10) %>% 
  mutate(Phi=round(Phi,4)) %>% 
  select(Shrub, Phi_round=Phi) %>% unique()

tol <- .0005

tdat <- pscores %>% left_join(., linedat, by="Shrub") %>% 
  filter(Phi <= (Phi_round + tol), Phi >= (Phi_round - tol))




###
# ITS 
###

#import
ps <- readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- data.frame(t(otu_table(ps) ))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")

#ID unique taxa, remove them
uasv <- (asvtab > 0)
uasv <- names(which(colSums(uasv)!=1))

asvtab <- asvtab %>% select(all_of(uasv))

#run program
isares.shrub <- multipatt(asvtab, meta$Shrub,
                          control=how(nperm = 1000), max.order = 2)
saveRDS(isares.shrub,"data/ITSIndSpec/res.shrub.ITS.rds")
#isares.shrub <- readRDS("data/ITSIndSpec/res.shrub.rds")

#interpret
spec <- data.frame(isares.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity")

sens <- data.frame(isares.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")

scores <- data.frame(isares.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")

gdat <- left_join(spec, sens, by=c("ASV","Shrub") ) %>% left_join(., scores, by=c("ASV","Shrub"))

tdat <- gdat %>% filter(gdat$Shrub %in% shrublevs )%>% 
  mutate(Shrub = factor(Shrub, levels = shrublevs  ))
g1b <- ggplot(tdat)+
  geom_function(fun=function(x){(0.75**2)/x}, color='red')+
  geom_point_rast(aes(x=specificity, y=sensitivity, fill=Shrub),alpha=0.8, shape=21, stroke=.25, size=3)+
  scale_y_continuous(expand=c(0,0), limits = c(0,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=shrubpal)+
  theme_bw()+
  xlab("Specificity")+
  ylab("Sensitivity")+
  ggtitle("Indicator Species Analysis (Microhabitat Single Level ITS)")
g1b
ggsave("figs/isa_shrub_simple.ITS.png", g1b)

#now with phi statistic
phires.shrub <- multipatt(asvtab, meta$Shrub, func = 'r.g',
                          control=how(nperm = 1000), max.order = 2)
saveRDS(phires.shrub,"data/ITSIndSpec/res.phi.shrub.ITS.rds")

pscores <- as.data.frame(phires.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))

g2b <- ggplot(pscores %>% filter(Shrub %in% c("Creosote","Mariola","Mesquite","Tarbush","Interspace" ) ) )+
  geom_jitter(aes(x = Shrub, y=Phi, fill=abs(Phi), alpha=abs(Phi) ), 
              height=0, shape=21,size=2)+
  scale_fill_gradientn(colors = c("black","darkred","darkorange","yellow"), name = "| Phi |")+
  scale_alpha_continuous(range =  c(0,0.9))+
  guides(alpha='none')+
  ggtitle("Phi Correlation (ITS OTUs across Shrubs)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
g2b
ggsave("figs/phi_shrub_simple_ITS.png", g2b)


#what about all shrubs vs interspace?
meta <- meta %>% mutate(ShrubCat = case_when(Shrub=="Interspace"~"Interspace",
                                             .default = "Shrub"))
#16S
ps <- readRDS("data/16Sdada2/phyloseq_dada2.OTU.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- data.frame(t(otu_table(ps)))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")

#ID unique taxa, remove them
uasv <- (asvtab > 0)
uasv <- names(which(colSums(uasv)!=1))

asvtab <- asvtab %>% select(all_of(uasv))

isa.shrub <- multipatt(asvtab, meta$Shrub,
                          control=how(nperm = 1000), min.order = 4)
saveRDS(isa.shrub,"data/16SIndSpec/res.shrub.interspace.16S.rds")
isares.shrub <- readRDS("data/16SIndSpec/res.shrub.16S.rds")

spec <- data.frame(isa.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity")
spec1 <- data.frame(isares.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity") %>% 
  filter(Shrub == "Interspace")

sens <- data.frame(isa.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")
sens1 <- data.frame(isares.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")%>% 
  filter(Shrub == "Interspace")

scores <- data.frame(isa.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")
scores1 <- data.frame(isares.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")%>% 
  filter(Shrub == "Interspace")

gdat <- left_join(spec, sens, by=c("ASV","Shrub") ) %>% left_join(., scores, by=c("ASV","Shrub"))
gdat1 <- left_join(spec1, sens1, by=c("ASV","Shrub") ) %>% left_join(., scores1, by=c("ASV","Shrub"))

tdat <- rbind(gdat,gdat1) %>% 
  filter(Shrub %in% c("Interspace", "Creosote.Mariola.Mesquite.Tarbush")) %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                              .default="Shrub"))
g1c <- ggplot(tdat)+
  geom_function(fun=function(x){(0.75**2)/x}, color='red')+
  geom_point_rast(aes(x=specificity, y=sensitivity, fill=Microhabitat),alpha=0.8, shape=21, stroke=.25, size=3)+
  scale_y_continuous(expand=c(0,0), limits = c(0,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=c(shrubpal[5],"#349D3B" )) +
  theme_bw()+
  xlab("Specificity")+
  ylab("Sensitivity")+
  ggtitle("Indicator Species Analysis (16S Core Microbiome)")
g1c
ggsave("figs/isa_shrub_v_interspace.16S.png", g1c)


phi.shrub <- multipatt(asvtab, meta$Shrub, func = 'r.g',
                       control=how(nperm = 1000), min.order = 4)
saveRDS(phi.shrub, "data/16SIndSpec/res.phi.shrub.interspace.16S.rds")
phires.shrub <- readRDS("data/16SIndSpec/res.phi.shrub.16S.rds")

pscores <- as.data.frame(phi.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                                       .default="Shrub"))
pscores1 <- as.data.frame(phires.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                                  .default="Shrub")) %>% filter(Shrub=="Interspace")

g2c <- ggplot(rbind(pscores,pscores1) ) +
  geom_jitter(aes(x = Microhabitat, y=Phi, fill=abs(Phi), alpha=abs(Phi) ), 
              height=0, shape=21,size=2)+
  scale_fill_gradientn(colors = c("black","darkred","darkorange","yellow"), name = "| Phi |")+
  scale_alpha_continuous(range =  c(0.2,0.9))+
  guides(alpha='none')+
  ggtitle("Phi Correlation (16S OTUs across Shrubs)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
g2c

ggsave("figs/phi_shrub_v_interspace.16S.png", g2c)

#ITS
ps <- readRDS("data/ITSdada2/phyloseq_dada2.OTU.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- data.frame(t(otu_table(ps)))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")

#ID unique taxa, remove them
uasv <- (asvtab > 0)
uasv <- names(which(colSums(uasv)!=1))

asvtab <- asvtab %>% select(all_of(uasv))

isa.shrub <- multipatt(asvtab, meta$Shrub,
                       control=how(nperm = 1000), min.order = 4)
saveRDS(isa.shrub, "data/ITSIndSpec/res.shrub.interspace.ITS.rds")
isares.shrub <- readRDS("data/ITSIndSpec/res.shrub.ITS.rds")
 
spec <- data.frame(isa.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity")
spec1 <- data.frame(isares.shrub$A) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "specificity") %>% 
  filter(Shrub == "Interspace")

sens <- data.frame(isa.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")
sens1 <- data.frame(isares.shrub$B) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sensitivity")%>% 
  filter(Shrub == "Interspace")

scores <- data.frame(isa.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")
scores1 <- data.frame(isares.shrub$str)%>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "sqrt.IndVal")%>% 
  filter(Shrub == "Interspace")

gdat <- left_join(spec, sens, by=c("ASV","Shrub") ) %>% left_join(., scores, by=c("ASV","Shrub"))
gdat1 <- left_join(spec1, sens1, by=c("ASV","Shrub") ) %>% left_join(., scores1, by=c("ASV","Shrub"))

tdat <- rbind(gdat,gdat1) %>% 
  filter(Shrub %in% c("Interspace", "Creosote.Mariola.Mesquite.Tarbush")) %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                                  .default="Shrub"))
g1d <- ggplot(tdat)+
  geom_function(fun=function(x){(0.75**2)/x}, color='red')+
  geom_point_rast(aes(x=specificity, y=sensitivity, fill=Microhabitat),alpha=0.8,shape=21, size=3, stroke=.25)+
  scale_y_continuous(expand=c(0,0), limits = c(0,1))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=c(shrubpal[5],"#349D3B" ))+
  theme_bw()+
  xlab("Specificity")+
  ylab("Sensitivity")+
  ggtitle("Indicator Species Analysis (ITS Core Microbiome)")
g1d
ggsave("figs/isa_shrub_v_interspace.ITS.png", g1d)


phi.shrub <- multipatt(asvtab, meta$Shrub, func = 'r.g',
                       control=how(nperm = 1000), min.order = 4)
saveRDS(phi.shrub, 'data/ITSIndSpec/res.phi.shrub.interspace.ITS.rds')
phires.shrub <- readRDS("data/ITSIndSpec/res.phi.shrub.ITS.rds")

pscores <- as.data.frame(phi.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                                  .default="Shrub"))
pscores1 <- as.data.frame(phires.shrub$str) %>% 
  rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = "Shrub", values_to = "Phi") %>% 
  mutate(Microhabitat = case_when(grepl("Interspace", Shrub)~"Interspace", 
                                  .default="Shrub")) %>% filter(Shrub=="Interspace")

g2d <- ggplot(rbind(pscores,pscores1) ) +
  geom_jitter(aes(x = Microhabitat, y=Phi, fill=abs(Phi), alpha=abs(Phi) ), 
              height=0, shape=21,size=2)+
  scale_fill_gradientn(colors = c("black","darkred","darkorange","yellow"), name = "| Phi |")+
  scale_alpha_continuous(range =  c(0.2,0.9))+
  guides(alpha='none')+
  ggtitle("Phi Correlation (ITS OTUs across Shrubs)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
g2d

ggsave("figs/phi_shrub_v_interspace.ITS.png", g2d)


library(patchwork)
png("figs/ISA_panel.png", width = 12, height = 15, units = 'in', res=500)
(g1a+g2a) /
  (g1b +g2b) /
  (g1c+g2c) /
  (g1d + g2d)
dev.off()




#write core microbiomes
out <- tdat %>% 
  filter(sqrt.IndVal>0.75,Microhabitat=="Shrub") %>% 
  left_join(., tax, by='ASV') %>% 
  mutate(across(Kingdom:Species, function(x){str_remove(x,'[a-z]__')})) %>% 
  select(Microhabitat, OTU=ASV, Sens=sensitivity, Spec=specificity, IndVal=sqrt.IndVal, Kingdom:Species)


write_csv(out, "data/ITS_shrub_core_microbiome.csv")
