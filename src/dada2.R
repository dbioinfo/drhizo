library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(colorspace)

setwd('~/WorkForaging/Academia/Nicole/drhizo')

#import metadata
reps <- read_csv('data/ShrubRhizo16S.env.processed.csv') %>% #use this file to remove duplicates
  select(Sample_ID=phinchID, Replicate) %>% 
  filter(Replicate=='Replicate') %>% 
  pull(Sample_ID)
meta <- read_csv('data/SraRunTable.csv') %>% 
  dplyr::rename(SampleID = "Run") %>% 
  filter(Sample_ID %in% reps, !(SampleID %in% c('SRR14874590','SRR14874591')) )%>% 
  select(SampleID, Site, Shrub) %>% 
  mutate(Shrub=factor(Shrub, levels=c("Interspace","Creosote","Mariola","Mesquite","Tarbush"))) %>% 
  as.data.frame()
meta <- rbind(meta, tibble(SampleID='SRR14874551',Site='AF2',Shrub='Creosote')) #fix oddball
rownames(meta) <- meta$SampleID


#list the forward and reverse reads in separate lists
fnFs <- sort(list.files('raw_data/fastq', pattern="_1.fastq", full.names = T))
fnRs <- sort(list.files('raw_data/fastq', pattern='_2.fastq', full.names=T))
snames <- gsub(pattern = '_1.fastq', x=list.files('raw_data/fastq', pattern='_1.fastq'),replacement = '')

#remove dups
fnFs <- fnFs[grepl(pattern = paste(meta$SampleID,collapse = '|'),x=fnFs)]
fnRs <- fnRs[grepl(pattern = paste(meta$SampleID,collapse = '|'),x=fnRs)]
snames <- snames[grepl(pattern = paste(meta$SampleID,collapse = '|'),x=snames)]

#check the quality of the reads, then adjust these settings
#plotQualityProfile(fnFs[1:3]); plotQualityProfile(fnRs[1:3])
trunc_f <- 250
trunc_r <- 160
left_trim <- 20
maxEE_f <- 2
maxEE_r <- 5

#filter and trim the data to ensure high quality reads only. As it stands, ~50% of reads are thrown out
filtFs <- file.path("raw_data/filt_fastq", paste0(snames, "_1.fq.gz"))
filtRs <- file.path("raw_data/filt_fastq", paste0(snames, "_2.fq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(trunc_f, trunc_r), trimLeft = left_trim, #rev reads terrible starting quality
                     maxN=0, maxEE=c(maxEE_f,maxEE_r), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
filt_fs <- list.files('raw_data/filt_fastq/', pattern='_1.fq.gz', full.names = T)
filt_rs <- list.files('raw_data/filt_fastq/', pattern='_2.fq.gz', full.names = T)
write.csv(out, 'data/dada2_reads_filt.csv', row.names = T)

#summarize filter stats 
out <- data.frame(out)
out$SampleID <- rownames(out)
out <- out %>% mutate(SampleID = str_remove(SampleID, "_[12]\\.fastq$"))
filtstat <- left_join(out, meta, by='SampleID') 
filtsum <- filtstat %>% 
  group_by(Shrub, Site) %>% 
  summarize(mean_reads_in=mean(reads.in), 
            mean_reads_filt=mean(reads.out), 
            sd_reads_in=sd(reads.in), 
            sd_reads_filt=sd(reads.out))
write_csv(filtsum, 'data/dada2_filt_stats.csv')
ggplot(filtstat)+ #reads in / out on y axis
  geom_boxplot(mapping=aes(x=Shrub, y=reads.out, fill=Shrub),outliers = F)+
  geom_point(mapping=aes(x=Shrub,y=reads.out))+
  facet_wrap(~Site)+
  theme_bw()+
  guides(fill='none')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = .5))
ggsave('figs/dada2/post_filt_reads.png')

#learn error rates (~10 min run time w/32 cores)
err_f <- learnErrors(filt_fs, multithread = T)
err_r <- learnErrors(filt_rs, multithread = T)

#run dada2
dadaFs <- dada(filt_fs, err_f, multithread = T)
dadaRs <- dada(filt_rs, err_r, multithread = T)

intout <- list(err_f=err_f, err_r=err_r, dadaFs=dadaFs, dadaRs=dadaRs)
saveRDS(intout, 'data/dada2_errs.rds')

#merge results together
mergers <- mergePairs(dadaFs, filt_fs, dadaRs, filt_rs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
#save
saveRDS(mergers,'data/merged_reads_dada2.rds')

#check some more summary stats
mergestat <- bind_rows(
  lapply(names(mergers), function(samplename) {
    dat <- mergers[[samplename]]
    tibble(
      SampleID     = gsub(pattern="_[12]\\.fq\\.gz$", replacement = '', samplename),
      n_unique_asvs = nrow(dat),
      reads.merged   = sum(dat$abundance),
      avg_overlap   = if (nrow(dat) > 0) mean(dat$nmatch) else NA_real_
    )
  })
)
mergestat <- left_join(mergestat, filtstat, by='SampleID') %>% na.omit()
ggplot(mergestat)+ #reads merged / unique asv on y axis
  geom_boxplot(mapping=aes(x=Shrub, y=n_unique_asvs, fill=Shrub),outliers = F)+
  geom_point(mapping=aes(x=Shrub,y=n_unique_asvs))+
  facet_wrap(~Site)+
  theme_bw()+
  guides(fill='none')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = .5))
ggsave('figs/dada2/unique_ASVs.png')
write_csv(mergestat,'data/dada2_mergestat.csv')

#filter chimeras (small amount present)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
1 - sum(seqtab.nochim)/sum(seqtab) #what percent of reads were chimera/bimera

#assign taxa from SILVA
tax <- assignTaxonomy(seqtab.nochim, "raw_data/silva_ref/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
tax <- addSpecies(tax, "raw_data/silva_ref/silva_v138.2_assignSpecies.fa.gz")

#create phyloseq obj
rownames(seqtab.nochim) <- gsub('_1.fq.gz','',rownames(seqtab.nochim))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(tax))

#fix ASV names, exact seqs are annoying
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#save phyloseq
otu <- data.frame(otu_table(ps))
otu$Sample <- rownames(otu)
otu <- otu %>% relocate(Sample)
write_csv(otu,'data/dada2_otu.csv')
saveRDS(ps, 'data/phyloseq_dada2.rds')

#graph phylum abundance
ab <- otu_table(ps) %>% t()
asvs<-rownames(ab)
ab <- as.data.frame(as.matrix(ab))
ab$OTU <- asvs
tax <- tax_table(ps) 
asv <- rownames(tax)
tax <- as.data.frame(as.matrix(tax))
tax$OTU <- asv

#arrange data by phylum
ab <- ab %>% 
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
ggplot(ab)+
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
ggsave('figs/dada2/phylum_abundance_plot.png', units = 'px', height=2000, width=5000)

#calculate alpha
alphas <- estimate_richness(ps, split = T, measures = c('Observed','Chao1','Simpson','Shannon'))
alphas$SampleID <- rownames(alphas)
write_csv(alphas, 'data/dada2_alphas.csv')

#graph alpha by shrub / facet site
gdat <- left_join(alphas, meta, by='SampleID')
ggplot(gdat)+
  geom_boxplot(aes(x=Shrub, y=Observed, fill=Shrub))+
  geom_jitter(aes(x=Shrub, y=Observed), height = 0, width = 0.05)+
  theme_bw()+
  facet_wrap(~Site)+
  guides(fill='none')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = .5))+
  ggtitle('Observed Entropy (Dada2)')
ggsave('figs/dada2/observed_alpha.png')




#output some summary stats
  #use below if starting from savepoint
  #mergers <- readRDS('data/merged_reads_dada2.rds')
  #intout <- readRDS('data/dada2_errs.rds')
  #list2env(intout, envir=.GlobalEnv)
getN <- function(x) sum(getUniques(x))
track <- data.frame(cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
track$SampleID <- rownames(track)
track <- track %>% mutate(SampleID=str_remove(SampleID, "_[12]\\.fq\\.gz$"))
track <- left_join(track, meta, by='SampleID')
track <- left_join(track, out %>% rename('reads.out'='reads.filtered'), by='SampleID')
track <- track %>% relocate(SampleID, Site, Shrub, reads.in, reads.filtered)
head(track)
write_csv(track, 'data/read_tracking.csv')
