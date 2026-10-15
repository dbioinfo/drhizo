# Desert Rhizosphere Microbiome

This repo contains several analyses completed for the Pietrasiask Lab pertaining to the rhizosphere microbiome of several shrub species in the mojave desert.

The repo is split into two main components -- `Data Processing` and `Data Analysis`. Whereas Data Processing comprises the transformation from raw reads into abundances, taxonomy assignment, and phylogenetic placement, Data Analysis comprises all the downstream analyses of these data products. 

## Data Processing Pipeline
```mermaid
graph TB;
    A[Raw Sequences: \nFASTQ] --> B[dada2: \nDenoising, ASV calling, \nFilter samples, reads & taxa];
    B --> C[dada2 (SILVA | UNITE): \nAssign Taxa];
    C --> D[qiime (SEPP): \nPhylogenetic Placement];
    D --> H[ape: \nPCoA];
    D --> G[indicspecies: \nIndicator Species Analysis];
    D --> D[vegan: \nAlpha diversity];
    D --> E[mecocomp (+SpiecEasi): \nBeta diversity];
    E --> F[iCAMP: \nCommunity Assembly];
    D --> I[faprotax | funguild : \nPathway Analysis];
```


## Data Analysis Pipeline
```mermaid
graph TB;
    A[Phyloseq object : \nASV table, Taxonomy, \nFunctional annotation, \nPhylo tree] --> B[`drhizo_16S_icamp_eda_shiny.Rmd`: \niCAMP Interactive Analysis];
    A --> C[`final_report_figs.R`  \nGenerate reports for: \nSample Distribution \nTaxonomic Abundance \nFunctional Annotation \nPCoA & Beta Div \nIndicator Species Analysis];
    A --> D['vegan_diversity.R' : \nAlpha Diversity Panel];
```

## Misc Scripts

The `asv_2_otu.R` script clusters ASVs determined by dada2 into OTUs with 95% similarity using the TreeLine function from DECIPHER. Interestingly, the liberal cutoff applied still left about half of 16S and ITS reads in a unique bin.

The scripts `gen_newick.sh` and `add_phylogeny.R` are a pair of programs that first add the ASVs/OTUs to a phylogenetic tree using qiime's SEPP algorithm, then adding that information back into the phyloseq object. 


 
