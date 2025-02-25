# Unsupervised-machine-learning-definition-of-the-microbial-rare-biosphere
(page under work)

Source code for the paper "Definition of the microbial rare biosphere through unsupervised machine learning" by Pascoal et al. (2024), in peer-review in Communications Biology journal.

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
