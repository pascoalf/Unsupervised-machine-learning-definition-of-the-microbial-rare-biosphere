# Packages
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(dplyr)
library(stringr)

#
selected_samples_mosj_2019 <- read.table("./data/selected_samples", header = TRUE)

#
samples_df_full <- data.frame(Sample = list.files("./fullLengthMOSJ2019/sequences")) %>% 
  mutate(Sample_id = str_extract(Sample, "M19_\\d+")) %>% 
  filter(Sample_id %in% selected_samples_mosj_2019$Sample)

#
seq_files_full_mosj <- paste0("./fullLengthMOSJ2019/sequences/", samples_df_full$Sample)

hist(nchar(getSequences("fullLengthMOSJ2019/sequences")), 100, 
     xlab = "Read length (bp)",
     main = "Full-length 16S rRNA gene raw reads size distribution")

# Primers, change as necessary
F27 <- "AGRGTTYGATYMTGGCTAG"
R1492 <- "RGYTACCTTGTTACGACTT"

# Make file for sequences with removed primers (removed in next step)
no_primers <- file.path("input", "noprimers", basename(seq_files_full_mosj))

# Remove primers
prim <- removePrimers(seq_files_full_mosj,
                      no_primers,
                      primer.fwd = F27,
                      primer.rev = dada2:::rc(R1492),
                      orient = TRUE)

# Make file path for filtered reads (filtered in next step)
filts_full_mosj <- file.path("./input", "noprimers", "filtered", basename(seq_files_full_mosj))

# Filter reads
# Note: Slow step (takes minutes).
track_full_mosj <- filterAndTrim(no_primers, ## verify this step
                                 filts_full_mosj,
                                 minQ = 3,
                                 minLen = 1000, 
                                 maxLen = 1600, 
                                 maxN = 0, 
                                 rm.phix = FALSE, 
                                 maxEE = 2, multithread = 10)

# Dereplicate reads
drp_full_mosj <- derepFastq(filts_full_mosj, verbose=TRUE)

# Learn error rates
# Note: very slow step (takes > 1h).
error_rates_full_mosj <- learnErrors(drp_full_mosj,
                                     errorEstimationFunction = PacBioErrfun,
                                     BAND_SIZE = 32,
                                     multithread = TRUE)

# Sanity check
#plotErrors(error_rates_full_mosj)

# Denoise reads
dd_seq_full_mosj <- dada(drp_full_mosj,
                         error_rates_full_mosj,
                         BAND_SIZE = 32,
                         multithread = TRUE) # Change to FALSE is you are using Windows OS.

# Make matrix with ASVs
seq_table_full_mosj <- makeSequenceTable(dd_seq_full_mosj)

# Remove chimeras
seq_tab_nochim_full_mosj <- removeBimeraDenovo(seq_table_full_mosj,
                                               minFoldParentOverAbundance = 3.5,
                                               multithread = TRUE)

# Track reads in each step
tracker_full_mosj <- cbind(ccs = prim[,1], 
                           primers = prim[,2],
                           filtered = track_full_mosj[,2],
                           denoised = lapply(dd_seq_full_mosj, function(x){sum(x$denoised)}),
                           nonchim =  apply(seq_tab_nochim_full_mosj,1,sum))

write.csv(tracker_full_mosj, "./data/tracker_full_mosj")

# Assign full taxonomy with Silva v138 database
ASV_taxonomy_silva_full_mosj <- assignTaxonomy(seq_tab_nochim_full_mosj,
                                               "~/Database/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                                               multithread = TRUE, 
                                               minBoot = 50)

# Assign full taxonomy with GTDB r202 database
ASV_taxonomy_gtdb_full_mosj <- assignTaxonomy(seq_tab_nochim_full_mosj, 
                                              "~/Database/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz",
                                              multithread = TRUE,
                                              minBoot = 50)
#
library(tidyr)

# Organize ASV abundance table
ASV_abundance_full_mosj <- seq_tab_nochim_full_mosj %>% t() %>% as.data.frame()
ASV_abundance_full_mosj$Sequence <- rownames(ASV_abundance_full_mosj)
rownames(ASV_abundance_full_mosj) <- NULL

# Silva taxonomy (full)
ASV_taxonomy_silva_full_mosj_df <- data.frame(ASV_taxonomy_silva_full_mosj)
ASV_taxonomy_silva_full_mosj_df$Sequence <- rownames(ASV_taxonomy_silva_full_mosj_df)
rownames(ASV_taxonomy_silva_full_mosj_df) <- NULL

ASV_abundance_taxonomy_full_mosj_df_silva <- ASV_abundance_full_mosj %>% 
  left_join(ASV_taxonomy_silva_full_mosj_df, by = "Sequence") %>% 
  mutate(Sequencing_technology = "PacBio",
         Database = "Silva") %>% 
  mutate(Species = ifelse(is.na(Species), NA, paste(Genus, Species)))

# GTDB taxonomy (full)
ASV_taxonomy_gtdb_full_mosj_df <- data.frame(ASV_taxonomy_gtdb_full_mosj)
ASV_taxonomy_gtdb_full_mosj_df$Sequence <- rownames(ASV_taxonomy_gtdb_full_mosj_df)
rownames(ASV_taxonomy_gtdb_full_mosj_df) <- NULL

ASV_abundance_taxonomy_full_mosj_df_gtdb <- ASV_abundance_full_mosj %>% 
  left_join(ASV_taxonomy_gtdb_full_mosj_df, by = "Sequence") %>% 
  mutate(Sequencing_technology = "PacBio",
         Database = "GTDB")

#
ASV_full_length_mosj <- ASV_abundance_taxonomy_full_mosj_df_silva %>% 
  rbind(ASV_abundance_taxonomy_full_mosj_df_gtdb)

#
mosj_2019_pacbio_tidy <- ASV_full_length_mosj %>% 
  #mutate(Species = ifelse(is.na(Species), NA, paste(Genus, Species))) %>% 
  pivot_longer(cols = contains("M19"), names_to = "Sample", values_to = "Abundance") %>% 
  pivot_longer(cols = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
               names_to = "Taxonomic_level", values_to = "Taxonomy")

write.csv(mosj_2019_pacbio_tidy, file = "./data/mosj_2019_pacbio_tidy")
