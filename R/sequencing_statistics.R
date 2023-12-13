## study statistics ##
## Arctic ocean dataset ##


# Helper functions

summarise_reads <- function(x){
  x %>% 
    group_by(Sample) %>% 
    filter(Abundance > 0) %>% 
    summarise(Reads = sum(Abundance))  %>% 
    summarise(meanReads = mean(Reads),
              sdReads = sd(Reads))  
}

summarise_asvs <- function(x){
  x %>% 
    group_by(Sample) %>% 
    filter(Abundance > 0) %>% 
    summarise(ASVs = specnumber(Abundance)) %>% 
    summarise(meanASVs = mean(ASVs),
              sdASVs= sd(ASVs))
}

## N-ICE dataset

# Initial high-quality reads and ASVs

# N-ICE ASVs
nice_asv %>% summarise_reads() 
nice_asv %>% summarise_asvs()

# N-ICE OTUs
nice_otu_raw %>% 
  prepare_tidy_data(sample_names = colnames(nice_otu_raw[,2:10])) %>% summarise_reads()
nice_otu_raw %>% 
  prepare_tidy_data(sample_names = colnames(nice_otu_raw[,2:10])) %>% summarise_asvs()

# N-ICE mOTUs
nice_mOTU_raw %>% 
  prepare_tidy_data(sample_names = colnames(nice_mOTU_raw[,2:10])) %>% summarise_reads()
nice_mOTU_raw %>% 
  prepare_tidy_data(sample_names = colnames(nice_mOTU_raw[,2:10])) %>% summarise_asvs()


# After pre-processing
nice_asv_clean %>% summarise_reads() 
nice_asv_clean %>% summarise_asvs()

# N-ICE OTUs
nice_otu %>% summarise_reads()
nice_otu %>% summarise_asvs()

# N-ICE mOTUs
nice_mOTU %>% summarise_reads()
nice_mOTU %>% summarise_asvs()

# After rarefaction
set.seed(123)
nice_dataset %>% 
  filter(Abundance > 1) %>% 
  mutate(Rarefaction = case_when(Type == "ASV" ~ asv_rarefaction,
                                 Type == "OTU" ~ otu_rarefaction,
                                 Type == "mOTU" ~ motu_rarefaction)) %>%
  group_by(Type, Sample, Rarefaction) %>% 
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(.x$Abundance, 
                                          sample = Rarefaction))))) %>% 
  unnest(c(data, Rarefied_reads)) %>% 
  rename(Rarefied_abundance = "V1") %>% 
  ungroup() %>% 
  group_by(Type, Sample) %>%
  filter(Rarefied_abundance > 0) %>% 
  summarise(ASVs = specnumber(Rarefied_abundance)) %>% 
  summarise(meanASVs = mean(ASVs),
            sdASVs= sd(ASVs))


## V4V5 vs full-length 16S
# Before pre-processing
# V4V5 reads
short_asv %>% 
  mutate(ID = paste0("ASV_", row_number(.))) %>% 
  pivot_longer(cols = contains("M19"), values_to = "Abundance", names_to = "Sample") %>% 
  group_by(Marker, Sample) %>% 
  summarise(HighQualityReads = sum(Abundance)) %>% 
  group_by(Marker) %>% 
  summarise(meanHighQualityReads = mean(HighQualityReads),
            sdHighQualityReads = sd(HighQualityReads))
# V4V5 asvs
short_asv %>% 
  mutate(ID = paste0("ASV_", row_number(.))) %>% 
  pivot_longer(cols = contains("M19"), values_to = "Abundance", names_to = "Sample") %>% 
  group_by(Marker, Sample) %>% 
  filter(Abundance> 0) %>% 
  summarise(ASVs = specnumber(Abundance)) %>% 
  summarise(meanASVs = mean(ASVs),
            sdASVs = sd(ASVs))

# full-length reads
ASV_abundance_taxonomy_full_mosj_df_silva %>% 
  mutate(ID = paste0("ASV_", row_number())) %>% 
  pivot_longer(cols = contains("M19"), values_to = "Abundance", names_to = "Sample") %>%
  select(ID, Sample, Marker, Abundance) %>% 
  mutate(Sample = str_remove(Sample, ".hifi_reads.fastq")) %>% 
  group_by(Marker, Sample) %>% 
  summarise(HighQualityReads = sum(Abundance)) %>% 
  group_by(Marker) %>% 
  summarise(meanHighQualityReads = mean(HighQualityReads),
            sdHighQualityReads = sd(HighQualityReads))
# full-length asvs
full_before_pre_processing %>% 
  group_by(Marker, Sample) %>% 
  summarise(ASVs = specnumber(Abundance)) %>% 
  group_by(Marker) %>% 
  summarise(meanASVs = mean(ASVs),
            sdASVs = sd(ASVs))

# After pre-processing
# both
# reads
short_and_full_asv %>% 
  group_by(Marker, Sample) %>% 
  summarise(HighQualityReads = sum(Abundance)) %>% 
  group_by(Marker) %>% 
  summarise(meanHighQualityReads = mean(HighQualityReads),
            sdHighQualityReads = sd(HighQualityReads))
## ASVs
short_vs_full_ulrb %>% 
  group_by(Marker, Sample) %>% 
  summarise(ASVs = specnumber(Abundance)) %>% 
  group_by(Marker) %>% 
  summarise(meanASVs = mean(ASVs),
            sdASVs = sd(ASVs))


## MOSJ 2016-2020 dataset
# ASVs after pre-processing 
all_years %>% 
  group_by(year, Sample) %>% 
  summarise(ASVs = specnumber(Abundance)) %>% 
  group_by(year) %>% 
  summarise(meanASVs = round(mean(ASVs),1),
            sdASVs = round(sd(ASVs),1)) 

# ASVs after rarefaction
all_years_rarefaction %>% 
  filter(Rarefied_abundance > 0) %>% 
  group_by(year, Sample) %>% 
  summarise(ASVs = specnumber(Rarefied_abundance)) %>% 
  group_by(year) %>% 
  summarise(meanASVs = round(mean(ASVs),1),
            sdASVs = round(sd(ASVs),1)) 


# total samples used (after all clean up)
all_years %>% 
  group_by(year) %>% 
  select(Sample) %>% 
  distinct() %>% 
  count()

# Initial reads
all_years_raw %>%
  group_by(year, Sample) %>% 
  summarise(Reads = sum(Abundance)) %>% 
  group_by(year) %>% 
  summarise(totalReads = round(mean(Reads),1),
            sdASVs = round(sd(Reads),1))

# Reads after pre-processing
all_years %>%
  group_by(year, Sample) %>% 
  summarise(Reads = sum(Abundance)) %>% 
  group_by(year) %>% 
  summarise(totalReads = round(mean(Reads),1),
            sdASVs = round(sd(Reads),1)) 

