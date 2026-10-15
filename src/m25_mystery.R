library(tidyverse)
library(phyloseq)
library(vegan)
library(stringdist)

#what is going on with M25?
#large number of unique taxa leading to high diversity measurements

#import data
setwd('~/WorkForaging/Academia/Nicole/drhizo')
ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')
otu <- as.data.frame(t(otu_table(ps)))
pres <- otu > 0


#number of unique asvs in M25 is 790
sid <- 'M25'
length(which(pres[,sid]))

#summary of number of unique asvs in every other sample. Max 374, less than half m25
summary(colSums(pres[,2:ncol(pres)]))

#how many asvs are unique to a sample? M25 has 309, the second highest is 138 and its a pretty normal distribution
asv_occurence <- rowSums(pres)
unique_asvs <- pres[asv_occurence == 1, , drop = FALSE]
unique_counts_per_sample <- colSums(unique_asvs)
summary(unique_counts_per_sample)
summary(unique_counts_per_sample[2:length(unique_counts_per_sample)])
unique_asv_list <- apply(unique_asvs, 2, function(x) rownames(unique_asvs)[x])

#which asvs are specific to m25?
unique_asv_m25 <- unique_asv_list[['M25']]

#are these ASVs super low abundance -- possible technical error
m25_abund <- otu[unique_asv_m25,'M25']
summary(m25_abund) #none of these are even single counts, some relatively high abundance

#compared to everything, how do these asvs fall into the distribution?
abnd <- rowSums(otu[, , drop=FALSE])
prev <- rowSums(otu>0)
uniq_label <- rownames(otu) %in% unique_asv_m25
tdf <- data.frame(abundance=abnd, prevalence=prev, m25_unique=uniq_label)

ggplot(tdf)+ #these asvs appear well within the distribution in this range...
  geom_point(aes(x=prevalence, y=abundance, color=m25_unique, size=m25_unique), alpha=.4)+
  scale_y_log10()+
  scale_size_manual(values = c(2,5))+
  theme_bw()

#maybe the sequences are really similar to another asvs? not particularly. 
seqs <- as.character(refseq(ps)) 
head(seqs[unique_asv_m25])
shared_asv_m25 <- (!(names(seqs) %in% unique_asv_m25) & pres[,'M25'])
dmat <- stringdistmatrix(
  seqs[unique_asv_m25],
  seqs[shared_asv_m25],
  method = "hamming"
)
hd_df <- data.frame(
  ASV = unique_asv_m25,
  min_hamming_to_shared = min_dist
)
summary(hd_df)
ggplot(hd_df)+ #the peak min hamming distance is around 8, likely not a technical error
  geom_histogram(aes(x=min_hamming_to_shared), binwidth = 1, fill='darkblue')+
  ggtitle("Number of differences between m25 unique ASVs and closest shared ASV")+
  xlab("Minimum Hamming Distance")+
  theme_bw()


#strange.. it seems these ASVs are genuine, so who are they?
#graph phylum abundance
ab <- otu_table(ps) %>% t()
asvs<-rownames(ab)
ab <- as.data.frame(as.matrix(ab))
ab$OTU <- asvs
tax <- tax_table(ps) 
asv <- rownames(tax)
tax <- as.data.frame(as.matrix(tax))
tax$OTU <- asv
meta <- data.frame(sample_data(ps))

#arrange data by phylum
ab <- ab %>% 
  filter(OTU %in% unique_asv_m25) %>% 
  pivot_longer(-OTU, names_to = "SRAid", values_to = "count") %>% 
  filter(count>0) %>% 
  left_join(., tax %>% select(OTU, Phylum), by="OTU") %>% 
  group_by(SRAid, Phylum) %>% 
  summarize(count=sum(count)) %>%
  group_by(SRAid) %>%  
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>% 
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) 
ab <- ab %>% mutate(Phylum=case_when(percent<1~"<1% Abundance",
                                     .default = Phylum)) %>% 
  ungroup() %>% 
  group_by(SRAid, Phylum) %>% 
  summarize(count=sum(count)) %>%
  ungroup() %>% 
  group_by(SRAid) %>% 
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>%   
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) %>% 
  left_join(meta %>% select(Site, Shrub, SRAid=SampleID), by='SRAid') 


palette <- qualitative_hcl(18, palette = "Dynamic") # or "Set3", "Classic Dark"
plevs <- unique(ab$Phylum)
palette <- setNames(palette, nm = plevs)

#add sampnames
tmp <- read_csv('data/SraRunTable.csv')  %>% select(SRAid = "Run", SampleName=`Sample Name`)
ab <- left_join(ab, tmp, by='SRAid')
#fix new guy
ab[is.na(ab$SampleName),'SampleName'] <- 'M25'
dim(ab)
ggplot(ab)+ # so it looks like these bacteria come from representative taxa as well. 
  geom_bar(aes(x=SampleName, y=percent, fill=Phylum), stat='identity', na.rm = T) +
  facet_wrap(~Site+Shrub, nrow=1, scales = "free_x",
             labeller = labeller(label_wrap_gen(width = 50)))+
  scale_fill_manual(values = palette)+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("Sample")+
  ylab("Relative Abundance")


#only thing I can think of is the fact that m25 had better read quality than any other sample, so more taxa are kept
plotQualityProfile(fl = 'raw_data/fastq/M25_1.fastq')
plotQualityProfile(fl = 'raw_data/fastq/M25_2.fastq')

#random normal sample, reverse reads are FAR worse, but there's twice as many...
plotQualityProfile(fl = 'raw_data/fastq/SRR14874600_1.fastq')
plotQualityProfile(fl = 'raw_data/fastq/SRR14874600_2.fastq')
#another
plotQualityProfile(fl = 'raw_data/fastq/SRR14874604_1.fastq')
plotQualityProfile(fl = 'raw_data/fastq/SRR14874604_2.fastq')
