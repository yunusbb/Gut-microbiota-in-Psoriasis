# Gut-microbiota-in-Psoriasis
This is a repository for scripts which were used to analyse human gut metagenomic data of psoriasis patients and healthy controls



0. 0_tripletones_removing.ipynb - dedicated to filtering psq object to remove singletones, doubletones, tripletones, which are used in subsequent analyses.
1. 1_diversity.ipynb - Computing alpha- and beta- diversity
2. 2_ordinates.ipynb - ordination for all the taxonomical levels
3. 2.5_PCA.ipynb - PCA (with zero inputation + clr transformation) on different taxonomic levels
4. 3_matching_pairs_PC.ipynb -  Matching pairs analysi
5. Wilcoxon_Sign_Rank_Test_v2.R - body of the function for wilcoxon paired analysis

### Power analysis:
1. 4_power_analysis_p1.ipynb - get bray pc and metadata
2. power_analysis_bray_pc_status.ipynb - power analysis + plot 

### Input
1. phyloseq-psor-control-metagenome.RData - psq object (not filtered, without normalization)


Here are described the steps of analysis. However, there are could by slight differences in packages and dataset versions resulted in slightly different numbers which do not affect the conclusions.
