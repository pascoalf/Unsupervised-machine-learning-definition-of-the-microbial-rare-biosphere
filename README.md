# Unsupervised-machine-learning-definition-of-the-microbial-rare-biosphere

Source code and source data for the paper "Definition of the microbial rare biosphere through unsupervised machine learning" by Pascoal et al. (2025), in peer-review in Communications Biology journal.

**Source data**

If you want to reproduce the figures and tables in the paper, please use the source data provided and run the scripts in the /R folder, following this order:

1 - R/prepareSession.R;

2 - R/phylogeneticUnitsComparisonNICE.R; (equivalent to section 1 Results)

3 - R/shortVsFullLength16S.R; (equivalent to section 2 Results)

4 - (equivalent to section 3 Results)

4.1 - R/prepare_mosj_data.R;

4.2 - R/sample_size_effect.R;

4.3 - R/species_number_effect.R;

4.4 - R/seq_power_effect.R;

4.5 - R/ulrb_vs_dataset_size_all.R

5 - (equivalent to section 4 Results)

5.1 - R/BCI_processing.R

5.2 - R/ulrb_vs_nonMicrobiome.R

6 - others

6.1 - R/rarefactionCurves.R

6.2 - R/sequencing_statistics.R

The DADA2 scripts to process 16S rRNA gene reads might take a long time to run and are available in the folder /DADA2. If you try to reproduce those steps, be careful with file paths, names, etc.

Raw FASTQ files are available in ENA, see Data Availability statement of the article.

We reproduced our code in a Linux (Fedora v38) machine, R version 4.3.2

**Source data**

All files are in .rds format, which is very compact and compatible with R. To open the file in an R session, use the function readRDS("name of file").

Short description of source data files:

- nice_ASVs_long_format.rds (ASV table in tidy/long format from N-ICE dataset)
- ASV_taxonomy_silva_full_mosj.rds (ASV table from PacBio sequencing (full-length 16S rRNA gene), 2019, of MOSJ. Taxonomic assignments based on Silva reference)
- seq_tab_nochim_full_mosj.rds (raw ASV table used to generate ASV_taxonomy_silva_full_mosj.rds)
- seqtab.nochim_mosj.rds (raw ASV table from V4V5 16S rRNA gene sequencing)
- short_mosj_taxa_silva.rds (ASV taxonomic assignments of N-ICE V4V5 16S rRNA gene sequencing data, using Silva reference)
- AO_2016_2020_silva.rds (ASV table, with ASVs classified using the Silva reference database. Includes metadata. Table in tidy/long format)
- octocoral_microbiome.rds (ASV table adapted from Keller-Costa et al., (2021))

**References:**

Pascoal, F., Branco, P., Torgo, L., Costa, R., & Magalhães, C. (2025). Definition of the microbial rare biosphere through unsupervised machine learning. Communications Biology (in Peer-Review).

Keller-Costa, T., Lago-Lestón, A., Saraiva, J. P., Toscan, R., Silva, S. G., Gonçalves, J., Cox, C. J., Kyrpides, N., Nunes da Rocha, U., & Costa, R. (2021). Metagenomic insights into the taxonomy, function, and dysbiosis of prokaryotic communities in octocorals. Microbiome, 9(1), 72. https://doi.org/https://doi.org/10.1186/s40168-021-01031-y



