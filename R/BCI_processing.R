## BCI data preparation
# get tsv file in https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
BCI_raw <- read.table("./source_data/FullMeasurementBCI.tsv", sep = "\t", header = TRUE)


# prepare data for count table
BCI_relevant_data <- BCI_raw %>% 
  select(SpeciesID, PlotCensusNumber, ExactDate, Status) %>%
  filter(Status == "alive") %>% 
  mutate(ExactDate = as_date(ExactDate),
         year = year(ExactDate)) %>% 
  filter(!is.na(year)) %>% 
  mutate(Sample = paste(PlotCensusNumber, year, sep = "_"))
  
# calculate counts
BCI_counts <- BCI_relevant_data %>%
  group_by(Sample) %>% 
  count(SpeciesID)

#
BCI_samples_report <- table(BCI_counts$Sample) %>% as.data.frame() %>% rename(Sample = Var1)

# Good samples
BCI_samples_to_keep <- BCI_samples_report %>% filter(Freq > 2) %>% pull(Sample)

# apply ulrb
BCI_ulrb <- BCI_counts %>% 
  filter(Sample %in% BCI_samples_to_keep) %>% 
  rename(Abundance = n) %>% 
  define_rb()

## General statistics
BCI_counts$Sample %>% unique() %>% length()
BCI_counts$SpeciesID %>% unique() %>% length()
##



