## DADA2 MOSJ combined dataset 2016-2020
init <- Sys.time()
## Combine all years
## reads denoising
library(dada2)
library(dplyr)
library(stringr)

## reads denoising
# 2016
#get path with sequences
seq_path_2016 <- "~/SequenceArchive/MOSJ/MOSJ2016/"

# get file names
fnFs_2016 <- paste(seq_path_2016, list.files(seq_path_2016, pattern = "R1.fastq"),sep="")
fnRs_2016 <-  paste(seq_path_2016, list.files(seq_path_2016, pattern = "R2.fastq"),sep="")

# Extract sample names
sample.names_2016 <- fnRs_2016 %>% 
  str_remove("_R2.fastq.bz2") %>% 
  str_remove("~/SequenceArchive/MOSJ/MOSJ2016/")

# Place filtered files in filtered/ subdirectory
filtFs_2016 <- file.path(seq_path_2016, "filtered", paste0(sample.names_2016, "_R1.fastq.bz2")) ## change me
filtRs_2016 <- file.path(seq_path_2016, "filtered", paste0(sample.names_2016, "_R2.fastq.bz2")) ## change me
names(filtFs_2016) <- sample.names_2016
names(filtRs_2016) <- sample.names_2016

## dada2 on mock community with optimal parameters
out_2016 <- filterAndTrim(fnFs_2016, filtFs_2016, fnRs_2016, filtRs_2016, 
                          truncLen=c(249,214),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE)

#Learn error rates
error_f_2016 <- learnErrors(filtFs_2016, multithread = TRUE) 
error_r_2016 <- learnErrors(filtRs_2016, multithread = TRUE) 

# get exact sequences
dada_Fs_2016 <- dada(filtFs_2016, err = error_f_2016, multithread = TRUE) 
dada_Rs_2016 <- dada(filtRs_2016, err = error_r_2016, multithread = TRUE) 

#
mergers_2016 <- mergePairs(dada_Fs_2016, filtFs_2016, dada_Rs_2016, filtRs_2016, verbose=TRUE)

# make sequence table
seqtab_2016 <- makeSequenceTable(mergers_2016)

# remove chimeras
seqtab.nochim_2016 <- removeBimeraDenovo(seqtab_2016,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2016 <- assignTaxonomy(seqtab.nochim_2016, 
                                      "~/Database/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

## 2017
#get path with sequences
seq_path_2017 <- "~/SequenceArchive/MOSJ/MOSJ2017/"

# get file names
fnFs_2017 <- paste(seq_path_2017, list.files(seq_path_2017, pattern = "R1.fastq"), sep="")
fnRs_2017 <-  paste(seq_path_2017, list.files(seq_path_2017, pattern = "R2.fastq"), sep="")

# Extract sample names
sample.names_2017 <- fnRs_2017 %>% 
  str_extract("S19-\\d+_\\d+")

# Place filtered files in filtered/ subdirectory
filtFs_2017 <- file.path(seq_path_2017, "filtered", paste0(sample.names_2017, "_R1.fastq")) ## change me
filtRs_2017 <- file.path(seq_path_2017, "filtered", paste0(sample.names_2017, "_R2.fastq")) ## change me
names(filtFs_2017) <- sample.names_2017
names(filtRs_2017) <- sample.names_2017

## dada2 on mock community with optimal parameters
out_2017 <- filterAndTrim(fnFs_2017, filtFs_2017, fnRs_2017, filtRs_2017, 
                          truncLen=c(249,214),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, 
                          multithread=TRUE,
                          trimLeft = c(19, 20)) 

#Learn error rates
error_f_2017 <- learnErrors(filtFs_2017, multithread = TRUE) 
error_r_2017 <- learnErrors(filtRs_2017, multithread = TRUE) 

# get exact sequences
dada_Fs_2017 <- dada(filtFs_2017, err = error_f_2017, multithread = TRUE) 
dada_Rs_2017 <- dada(filtRs_2017, err = error_r_2017, multithread = TRUE) 

#
mergers_2017 <- mergePairs(dada_Fs_2017, filtFs_2017, dada_Rs_2017, filtRs_2017, verbose=TRUE)

# make sequence table
seqtab_2017 <- makeSequenceTable(mergers_2017)

# remove chimeras
seqtab.nochim_2017 <- removeBimeraDenovo(seqtab_2017,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2017 <- assignTaxonomy(seqtab.nochim_2017, 
                                      "~/Database/silva_nr99_v138.1_train_set.fa.gz")

## 2018
#get path with sequences
seq_path_2018 <- "~/SequenceArchive/MOSJ/MOSJ2018/"

# get file names
fnFs_2018 <- paste(seq_path_2018, list.files(seq_path_2018, pattern = "R1.fastq"), sep="")
fnRs_2018 <-  paste(seq_path_2018, list.files(seq_path_2018, pattern = "R2.fastq"), sep="")

# Extract sample names
sample.names_2018 <- fnRs_2018 %>% 
  str_remove("_R2.fastq.bz2") %>% 
  str_remove("~/SequenceArchive/MOSJ/MOSJ2018/")

# Place filtered files in filtered/ subdirectory
filtFs_2018 <- file.path(seq_path_2018, "filtered", paste0(sample.names_2018, "_R1.fastq.bz2")) ## change me
filtRs_2018 <- file.path(seq_path_2018, "filtered", paste0(sample.names_2018, "_R2.fastq.bz2")) ## change me
names(filtFs_2018) <- sample.names_2018
names(filtRs_2018) <- sample.names_2018

## dada2 on mock community with optimal parameters
out_2018 <- filterAndTrim(fnFs_2018, filtFs_2018, fnRs_2018, filtRs_2018, 
                          truncLen=c(249,214),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE, 
                          trimLeft = c(19, 20)) 

#Learn error rates
error_f_2018 <- learnErrors(filtFs_2018, multithread = TRUE) 
error_r_2018 <- learnErrors(filtRs_2018, multithread = TRUE) 

# get exact sequences
dada_Fs_2018 <- dada(filtFs_2018, err = error_f_2018, multithread = TRUE) 
dada_Rs_2018 <- dada(filtRs_2018, err = error_r_2018, multithread = TRUE) 

#
mergers_2018 <- mergePairs(dada_Fs_2018, filtFs_2018, dada_Rs_2018, filtRs_2018, verbose=TRUE)

# make sequence table
seqtab_2018 <- makeSequenceTable(mergers_2018)

# remove chimeras
seqtab.nochim_2018 <- removeBimeraDenovo(seqtab_2018,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2018 <- assignTaxonomy(seqtab.nochim_2018,  "~/Database/silva_nr99_v138.1_train_set.fa.gz")

## 2019
#get path with sequences
seq_path_2019 <- "~/SequenceArchive/MOSJ/MOSJ2019/Illumina/"

# get file names
fnFs_2019 <- paste(seq_path_2019, list.files(seq_path_2019, pattern = "R1_001.fastq"), sep="")
fnRs_2019 <- paste(seq_path_2019, list.files(seq_path_2019, pattern = "R2_001.fastq"), sep="")

# Extract sample names
sample.names_2019 <- list.files(seq_path_2019, pattern = "R1") %>% 
  str_remove("_L001_R1_001.fastq") 

# Place filtered files in filtered/ subdirectory
filtFs_2019 <- file.path(seq_path_2019, "filtered", list.files(seq_path_2019, pattern = "R1")) ## change me
filtRs_2019 <- file.path(seq_path_2019, "filtered", list.files(seq_path_2019, pattern = "R2")) ## change me
names(filtFs_2019) <- sample.names_2019
names(filtRs_2019) <- sample.names_2019

## dada2 on mock community with optimal parameters
out_2019 <- filterAndTrim(fnFs_2019, filtFs_2019, fnRs_2019, filtRs_2019, 
                          truncLen = c(249,214),
                          maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                          compress = TRUE, multithread = TRUE,
                          trimLeft = c(19, 20)) 

#Learn error rates
error_f_2019 <- learnErrors(filtFs_2019, multithread = TRUE) 
error_r_2019 <- learnErrors(filtRs_2019, multithread = TRUE) 

# get exact sequences
dada_Fs_2019 <- dada(filtFs_2019, err = error_f_2019, multithread = TRUE) 
dada_Rs_2019 <- dada(filtRs_2019, err = error_r_2019, multithread = TRUE) 

#
mergers_2019 <- mergePairs(dada_Fs_2019, filtFs_2019, dada_Rs_2019, filtRs_2019, verbose=TRUE)

# make sequence table
seqtab_2019 <- makeSequenceTable(mergers_2019)

# remove chimeras
seqtab.nochim_2019 <- removeBimeraDenovo(seqtab_2019,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2019 <- assignTaxonomy(seqtab.nochim_2019, 
                                      "~/Database/silva_nr99_v138.1_train_set.fa.gz",
                                      multithread = TRUE)

## 2020
#get path with sequences
seq_path_2020 <- "~/SequenceArchive/MOSJ/MOSJ2020/"

# get file names
fnFs_2020 <- paste(seq_path_2020, list.files(seq_path_2020, pattern = "R1_001.fastq"), sep="")
fnRs_2020 <- paste(seq_path_2020, list.files(seq_path_2020, pattern = "R2_001.fastq"), sep="")

# Extract sample names
sample.names_2020 <- list.files(seq_path_2020, pattern = "R1") %>% 
  str_remove("_L001_R1_001.fastq") 

# Place filtered files in filtered/ subdirectory
filtFs_2020 <- file.path(seq_path_2020, "filtered", list.files(seq_path_2020, pattern = "R1")) ## change me
filtRs_2020 <- file.path(seq_path_2020, "filtered", list.files(seq_path_2020, pattern = "R2")) ## change me
names(filtFs_2020) <- sample.names_2020
names(filtRs_2020) <- sample.names_2020


## dada2 on mock community with optimal parameters
out_2020 <- filterAndTrim(fnFs_2020, filtFs_2020, fnRs_2020, filtRs_2020, 
                          truncLen=c(249,214),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE, 
                          trimLeft = c(19, 20)) 

#Learn error rates
error_f_2020 <- learnErrors(filtFs_2020, multithread = TRUE) 
error_r_2020 <- learnErrors(filtRs_2020, multithread = TRUE) 

# get exact sequences
dada_Fs_2020 <- dada(filtFs_2020, err = error_f_2020, multithread = TRUE) 
dada_Rs_2020 <- dada(filtRs_2020, err = error_r_2020, multithread = TRUE) 

#
mergers_2020 <- mergePairs(dada_Fs_2020, filtFs_2020, dada_Rs_2020, filtRs_2020, verbose=TRUE)

# make sequence table
seqtab_2020 <- makeSequenceTable(mergers_2020)

# remove chimeras
seqtab.nochim_2020 <- removeBimeraDenovo(seqtab_2020,
                                         method = "consensus",
                                         multithread = TRUE,
                                         verbose = TRUE)

# Taxonomic classification down to genus level (naive Bayesian classifier)
taxonomy_silva_2020 <- assignTaxonomy(seqtab.nochim_2020, 
                                      "~/Database/silva_nr99_v138.1_train_set.fa.gz",
                                      multithread = "TRUE")

