# Unsupervised-machine-learning-definition-of-the-microbial-rare-biosphere
Source code for the paper "Unsupervised machine learning definition of the microbial rare biosphere" by Pascoal et al. (2023, in preparation for submission).

If you want to reproduce the figures and tables in the paper, please use the source data provided and run the scripts in the /R folder, following this order:

1 - prepareSession.R (will load R packages, etc);
2 - phylogeneticUnitsComparisonNICE.R (equivalent to section 1)
3 - shortVsFullLength16S.R (equivalent to section 2)
4 - ecologyArcticOcean2016To2020.R (equivalent to section 3)
5 - ulrb_vs_FuzzyQ.R (equivalent to section 4)
6 - rarefactionCurves.R
7 - sequencing_statistics.R

The DADA2 scripts might take a long time and are available in the folder /DADA2. If you try to reproduce those steps, be careful with file paths, names, etc.

Raw FASTQ files will be made available in ENA repository.
