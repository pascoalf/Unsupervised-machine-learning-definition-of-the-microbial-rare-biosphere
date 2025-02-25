# mosj 2016-2020 dataset
all_years_raw <- readRDS("./source_data/AO_2016_2020_silva.rds")

## Clean taxonomy
all_years <- all_years_raw %>% 
  filter(!is.na(Kingdom),
         !is.na(Phylum)) %>%
  mutate(Phylum = ifelse(Phylum == "Cyanobacteria/Chloroplast", 
                         "Cyanobacteria", Phylum)) %>% 
  filter(Order != "Chloroplast")

## Samples that pass rarefaction
samples_to_keep <- 
  all_years %>% 
  group_by(Sample) %>% 
  summarise(Total_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  filter(Total_reads >= 10000) %>% 
  pull(Sample)

# will give warning, this is not a problem
set.seed(123); all_years_rarefaction <- 
  all_years %>% 
  filter(Sample %in% samples_to_keep) %>%  ## filter samples with more than 10 000 reads
  filter(Abundance > 1) %>% 
  group_by(Sample) %>%
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(
                                    .x$Abundance, 
                                    sample = 10000))))) %>% 
  unnest(c(data,Rarefied_reads)) %>% 
  rename(Rarefied_abundance = "V1")

# alternative normalization
all_years_rclr <- all_years %>% 
  filter(Abundance > 1) %>% 
  group_by(Sample) %>% 
  mutate(sample_mean = mean(Abundance),
         clr = log(Abundance/sample_mean))


## Rare Biosphere ##
all_years_ulrb <- all_years_rarefaction %>%  
  mutate(Abundance = Rarefied_abundance, 
         Rarefied_abundance = NULL,
         Rarefaction = "Yes") %>% 
  define_rb()