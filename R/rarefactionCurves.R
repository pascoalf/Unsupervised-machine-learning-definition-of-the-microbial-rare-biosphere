## Rarefaction curves ##

### N-ICE ###
nice_dataset_wide <- nice_dataset %>% 
  pivot_wider(names_from = ID, 
              values_from = Abundance,
              values_fn = list(count = list)) %>% 
  unchop(everything()) %>% 
  distinct()

#
make_rarefaction_phylo_units <- function(unit = "ASV", ...){
  x <- nice_dataset_wide %>% filter(Type == unit)
  sample_names_rarefaction <-x$Sample 
  rownames(x) <- x$Sample
  x$Sample <- NULL
  x[is.na(x)] <- 0 
  x <- x[, -c(1,2)]
  #
  rownames(x) <- sample_names_rarefaction
  rarecurve(x, ...)
}
par(mfrow = c(1,3))
# ASVs
make_rarefaction_phylo_units(step = 100, 
                             lwd = 2, 
                             col = qualitative_colors[1],
                             main = "Rarefaction curve for ASVs (N-ICE, 2015)",
                             ylab = "Number of ASVs",
                             xlab = "Number of reads")
# OTUs
make_rarefaction_phylo_units(unit = "OTU", 
                             step = 100, 
                             lwd = 2, 
                             col = qualitative_colors[2],
                             main = "Rarefaction curve for OTUs (N-ICE, 2015)",
                             ylab = "Number of OTUs",
                             xlab = "Number of reads")
# mOTU
make_rarefaction_phylo_units(unit = "mOTU", 
                             step = 10, 
                             lwd = 2, 
                             col = qualitative_colors[3],
                             main = "Rarefaction curve for mOTUs (N-ICE, 2015)",
                             ylab = "Number of mOTUs",
                             xlab = "Number of reads")

##

### MOSJ 2019 ###
### N-ICE ###
short_and_full_asv_wide <- short_and_full_asv %>% 
  pivot_wider(names_from = ID, 
              values_from = Abundance,
              values_fn = list(count = list)) %>% 
  unchop(everything()) %>% 
  distinct()

#
make_rarefaction_short_vs_full <- function(marker = "V4-V5", ...){
  x <- short_and_full_asv_wide %>% filter(str_detect(Marker, marker))
  sample_names <- x$Sample
  x[is.na(x)] <- 0 
  x$Sample <- NULL
  x <- x[, -c(1,2)]
  x <- as.matrix(x)
  rownames(x) <- sample_names
  #
  rarecurve(x, ...)
}
par(mfrow = c(1,2))
# V4V5
make_rarefaction_short_vs_full(step = 1000, 
                               lwd = 2.5, 
                               col = qualitative_colors[1],
                               main = "V4V5 16S rRNA (MOSJ2019)",
                               xlab = "Number of reads",
                               ylab = "Number of ASVs")
# full-legth
make_rarefaction_short_vs_full(marker = "Full-length", 
                               step = 100, 
                               lwd = 2.5, 
                               col = qualitative_colors[2],
                               main = "Full-length 16S rRNA (MOSJ2019)",
                               xlab = "Number of reads",
                               ylab = "Number of ASVs")

## Arctic ocean combined dataset

all_years_wide <- all_years %>% 
  ungroup() %>% 
  mutate(Sample_unique = ifelse(year != 2015, paste(Sample), Sample)) %>% 
  select(Abundance, Sequence, Sample_unique, year) %>% 
  pivot_wider(names_from = Sequence, 
              values_from = Abundance,
              values_fn = list(count = list)) %>% 
  unchop(everything()) %>% 
  distinct()

#
make_rarefaction_AO <- function(y = 2020, ...){
  x <- all_years_wide %>% 
    filter(year == y)
  sample_names <- x$Sample_unique
  x[is.na(x)] <- 0 
  x$Sample <- NULL
  x <- x[, -c(1,2)]
  x <- as.matrix(x)
  rownames(x) <- sample_names
  #
  rarecurve(x, ...)
}


par(mfrow = c(2,3))
# 2015
#make_rarefaction_AO(y = 2015, step = 500, lwd = 2, col = qualitative_colors[1],
#                   main = "2015", xlab = "Number of reads", ylab = "Number of ASVs")
# 2016
make_rarefaction_AO(y = 2016, step = 500, lwd = 2, col = qualitative_colors[2],
                    main = "2016", xlab = "Number of reads", ylab = "Number of ASVs")
# 2017
make_rarefaction_AO(y = 2017, step = 500, lwd = 2, col = qualitative_colors[3],
                    main = "2017", xlab = "Number of reads", ylab = "Number of ASVs")
# 2018
make_rarefaction_AO(y = 2018, step = 500, lwd = 2, col = qualitative_colors[4],
                    main = "2018", xlab = "Number of reads", ylab = "Number of ASVs")
# 2019
make_rarefaction_AO(y = 2019, step = 500, lwd = 2, col = qualitative_colors[5],
                    main = "2019", xlab = "Number of reads", ylab = "Number of ASVs")
# 2020
make_rarefaction_AO(y = 2020, step = 500, lwd = 2, col = qualitative_colors[6],
                    main = "2020", xlab = "Number of reads", ylab = "Number of ASVs")

