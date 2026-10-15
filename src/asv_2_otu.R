library(tidyverse)
library(phyloseq)
library(dada2)
library(DECIPHER)

setwd('~/WorkForaging/Academia/Nicole/drhizo/')

get_conseq = function(x){
  if (length(x) > 1){
    aln = AlignSeqs(DNAStringSet(x))
  } else {
    aln = DNAStringSet(x)
  }
  con = ConsensusSequence(aln)
  return(as.character(con))
}

###
# 16S 
###

#import
ps <- readRDS("data/16Sdada2/phyloseq_dada2.tree.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- data.frame(otu_table(ps))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")
seqtab <- data.frame(refseq(ps)) %>% rename("refseq.ps."="Seq") %>% rownames_to_column("ASV")

seqs <- tax %>% 
  filter(!is.na(Phylum),
         Kingdom %in% c('Bacteria','Archaea')) %>% 
  left_join(., seqtab, by="ASV")

aln <- AlignSeqs(DNAStringSet(seqs$Seq))
dmat <- DistanceMatrix(aln, type = 'dist')
clsts = TreeLine(myDistMatrix = dmat, cutoff = 0.05, 
                 type="clusters", method = 'UPGMA', 
                 processors = 4)
clsts$seq <- seqs$Seq
clsts$ASV <- seqs$ASV

#prevalence of each asv
uasv <- (asvtab > 0)
gdat <- data.frame(ASV=colnames(uasv), Prev=colSums(uasv)) 
ggplot(gdat)+
  geom_density(aes(x=Prev), fill="#00FF0040")+theme_bw()+
  scale_x_continuous(expand = c(0,0))+  scale_y_continuous(expand = c(0,0))


conseq <- clsts %>%
  group_by(cluster) %>% 
  summarise(consensus = get_conseq(seq))
clstax = assignTaxonomy(conseq$consensus, 
                        refFasta = "raw_data/silva_ref/silva_nr99_v138.2_toGenus_trainset.fa.gz",
                          tryRC = TRUE, multithread = 10)
sclstax <-  assignTaxonomy(conseq$consensus, 
                           refFasta = "raw_data/silva_ref/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
                           tryRC = TRUE, multithread = 10)
                          
asvdf <- t(asvtab) %>% data.frame() %>% 
  mutate(ASV = rownames(.)) %>% 
  select(ASV, everything())

asvdf <- asvdf %>% 
  left_join(clsts, by = 'ASV')%>% 
  select(cluster, everything())
clstdf <- asvdf %>% select(-ASV, -seq) %>% 
  group_by(cluster) %>% 
  summarise_all(sum) %>% ungroup() %>% 
  filter(!is.na(cluster)) %>% data.frame()
dim(clstdf)

clstdf <- clstdf %>% mutate(OTU = paste0("OTU", cluster) ) %>% column_to_rownames("OTU")
sclstax <- data.frame(sclstax) %>% rownames_to_column("consensus") %>% left_join(., conseq, by="consensus") %>% 
  mutate(OTU = paste0("OTU", cluster)) %>% column_to_rownames("OTU")

clseq <- conseq %>% mutate(OTU = paste0("OTU", cluster)) %>% column_to_rownames("OTU") %>% select(-cluster)
rseq <- DNAStringSet(clseq$consensus)  
names(rseq) <- rownames(clseq)

psc <- phyloseq(otu_table = otu_table( clstdf %>% select(-cluster), taxa_are_rows = T ),
                tax_table = tax_table(as.matrix(sclstax %>% select(-consensus, -cluster))) )
psc <- merge_phyloseq(psc, rseq)


meta <- read_csv('data/SraRunTable.csv') %>% 
  select(SRA=Run, Sample=Sample_ID,Shrub,Site) %>% 
  column_to_rownames("SRA") %>% 
  sample_data()
psc <- merge_phyloseq(psc, meta)

saveRDS(psc, "data/16Sdada2/phyloseq_dada2.OTU.rds")







###
# ITS 
###

#import
ps <- readRDS("data/ITSdada2/phyloseq_dada2.tree.rds")
meta <- data.frame(sample_data(ps)) %>% 
  mutate(Shrub = factor(Shrub, levels = c("Creosote","Mariola","Mesquite","Tarbush","Interspace" )  ))
asvtab <- data.frame(otu_table(ps))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")
seqtab <- data.frame(refseq(ps)) %>% rename("refseq.ps."="Seq") %>% rownames_to_column("ASV")

seqs <- tax %>% 
  filter(!is.na(Phylum),
         Kingdom %in% c('k__Fungi')) %>% 
  left_join(., seqtab, by="ASV")

aln <- AlignSeqs(DNAStringSet(seqs$Seq))
dmat <- DistanceMatrix(aln, type = 'dist')
clsts = TreeLine(myDistMatrix = dmat, cutoff = 0.05, 
                 type="clusters", method = 'UPGMA', 
                 processors = 8)
clsts$seq <- seqs$Seq
clsts$ASV <- seqs$ASV

#prevalence of each asv
uasv <- (asvtab > 0)
gdat <- data.frame(ASV=colnames(uasv), Prev=colSums(uasv)) 
ggplot(gdat)+
  geom_density(aes(x=Prev), fill="#00FF0040")+theme_bw()+
  scale_x_continuous(expand = c(0,0))+  scale_y_continuous(expand = c(0,0))


conseq <- clsts %>%
  group_by(cluster) %>% 
  summarise(consensus = get_conseq(seq))
clstax = assignTaxonomy(conseq$consensus, 
                        refFasta = "raw_data/unite_ref/sh_general_release_dynamic_19.02.2025.fasta",
                        tryRC = TRUE, multithread = 10)

asvdf <- t(asvtab) %>% data.frame() %>% 
  mutate(ASV = rownames(.)) %>% 
  select(ASV, everything())

asvdf <- asvdf %>% 
  left_join(clsts, by = 'ASV')%>% 
  select(cluster, everything())
clstdf <- asvdf %>% select(-ASV, -seq) %>% 
  group_by(cluster) %>% 
  summarise_all(sum) %>% ungroup() %>% 
  filter(!is.na(cluster)) %>% data.frame()
dim(clstdf)

clstdf <- clstdf %>% mutate(OTU = paste0("OTU", cluster) ) %>% column_to_rownames("OTU")
clstax <- data.frame(clstax) %>% rownames_to_column("consensus") %>% left_join(., conseq, by="consensus") %>% 
  mutate(OTU = paste0("OTU", cluster)) %>% column_to_rownames("OTU")

clseq <- conseq %>% mutate(OTU = paste0("OTU", cluster)) %>% column_to_rownames("OTU") %>% select(-cluster)
rseq <- DNAStringSet(clseq$consensus)  
names(rseq) <- rownames(clseq)

psc <- phyloseq(otu_table = otu_table( clstdf %>% select(-cluster), taxa_are_rows = T ),
                tax_table = tax_table(as.matrix(clstax %>% select(-consensus, -cluster))) )
psc <- merge_phyloseq(psc, rseq)

meta <- read_csv('data/SraRunTable.csv') %>% 
  select(SRA=Run, Sample=Sample_ID,Shrub,Site) %>% 
  filter(Sample %in% sample_names(psc)) %>% 
  column_to_rownames("Sample") %>% 
  sample_data()
psc <- merge_phyloseq(psc, meta)


saveRDS(psc, "data/ITSdada2/phyloseq_dada2.OTU.rds")

