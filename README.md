

# Desert Rhizosphere Microbiome

This repo contains several analyses completed for the Pietrasiask Lab pertaining to the rhizosphere microbiome of several shrub species in the mojave desert.

The repo is split into two main components -- `Data Processing` and `Data Analysis`. Whereas Data Processing comprises the transformation from raw reads into abundances, taxonomy assignment, and phylogenetic placement, Data Analysis comprises all the downstream analyses of these data products. 

## Pipeline Overview
```mermaid
---
config:
  theme: 'neo-dark'
---
graph TB;
subgraph Z[Data Processing];
    A((Raw Sequence Reads)) --> B[dada2: \nDenoising, ASV calling, \nFilter samples, reads & taxa];
    B --> C[dada2 SILVA & UNITE: \nAssign Taxa];
    C --> K[qiime SEPP: \nPhylogenetic Placement];
    K --> D((Complete Phyloseq Object));
end;
subgraph Y[Data Analysis];
    D((Complete Phyloseq Object)) --> H[ape: \nPCoA];
    D --> G[indicspecies: \nIndicator Species Analysis];
    D --> J[vegan: \nAlpha diversity];
    D --> E[mecocomp +SpiecEasi: \nBeta diversity];
    E --> F[iCAMP: \nCommunity Assembly];  
    D --> I[faprotax & funguild : \nPathway Analysis];
    G --> L[final_report_figs.R]; 
    J --> L;
    E --> L;
    F --> L;
    I --> L;
    H --> L;
    F --> N[`drhizo_16S_icamp_eda_shiny.Rmd`: \niCAMP Interactive Analysis];
    L --> M((final_report.pdf));
end;
```


## Misc Scripts

The `asv_2_otu.R` script clusters ASVs determined by dada2 into OTUs with 95% similarity using the TreeLine function from DECIPHER. Interestingly, the liberal cutoff applied still left about half of 16S and ITS reads in a unique bin.

The scripts `gen_newick.sh` and `add_phylogeny.R` are a pair of programs that first add the ASVs/OTUs to a phylogenetic tree using qiime's SEPP algorithm, then adding that information back into the phyloseq object. 


 
