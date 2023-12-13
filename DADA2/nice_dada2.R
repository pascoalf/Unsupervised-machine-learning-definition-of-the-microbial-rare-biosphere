library(dada2)
library(dplyr)
library(stringr)

## 2015
#get path with sequences
seq_path_2015 <- "~/SequenceArchive/N-ICE/16S/"

# get file names
fnFs_2015 <- paste(seq_path_2015, list.files(seq_path_2015, pattern = "_1.fastq.gz"), sep="")
fnRs_2015 <-  paste(seq_path_2015, list.files(seq_path_2015, pattern = "_2.fastq.gz"), sep="")

# Extract sample names
sample.names_2015 <- list.files(seq_path_2015, pattern = "_1.fastq.gz") %>% 
  str_remove("_2.fastq.gz")

# Place filtered files in filtered/ subdirectory
filtFs_2015 <- file.path(seq_path_2015, "filtered", paste0(sample.names_2015, "_1.fastq.gz")) ## change me
filtRs_2015 <- file.path(seq_path_2015, "filtered", paste0(sample.names_2015, "_2.fastq.gz")) ## change me
names(filtFs_2015) <- sample.names_2015
names(filtRs_2015) <- sample.names_2015

## dada2 on mock community with optimal parameters
out_2015 <- filterAndTrim(fnFs_2015, filtFs_2015, fnRs_2015, filtRs_2015, 
                          truncLen=c(249,214),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE)

#Learn error rates
error_f_2015 <- learnErrors(filtFs_2015, multithread = TRUE) 
error_r_2015 <- learnErrors(filtRs_2015, multithread = TRUE) 

# get exact sequences
dada_Fs_2015 <- dada(filtFs_2015, err = error_f_2015, multithread = TRUE) 
dada_Rs_2015 <- dada(filtRs_2015, err = error_r_2015, multithread = TRUE) 

#
mergers_2015 <- mergePairs(dada_Fs_2015, filtFs_2015, dada_Rs_2015, filtRs_2015, verbose=TRUE)

# make sequence table
seqtab_2015 <- makeSequenceTable(mergers_2015)

# remove chimeras
seqtab.nochim_2015 <- removeBimeraDenovo(seqtab_2015,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2015 <- assignTaxonomy(seqtab.nochim_2015, "~/Database/silva_nr99_v138.1_train_set.fa.gz")
taxonomy_silva_2015_df <- taxonomy_silva_2015 %>% 
  as.data.frame() %>% 
  mutate(Sequence = rownames(.))

## don't run
#save(list = ls(), file = "./intermediate_files/int_files_nice")
#load("./intermediate_files/int_files_nice")

nice_ASV_table <- seqtab.nochim_2015 %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Sequence = rownames(.)) %>% 
  left_join(taxonomy_silva_2015_df, by = "Sequence") %>% 
  mutate(ASV = paste0("ASV_", rownames(.)))

rownames(nice_ASV_table) <- NULL

nice_asvs_tidy <- prepare_tidy_data(nice_ASV_table, samples_in = "cols", sample_names = rownames(seqtab.nochim_2015)) %>% 
  mutate(Sample = str_remove(Sample, "_1.fastq.gz"))

saveRDS(nice_asvs_tidy, file = "./source_data/nice_ASVs_long_format.rds")
