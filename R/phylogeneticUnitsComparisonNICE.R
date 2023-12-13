## N-ICE dataset
# ulrb against phylogenetic units

## ASVs
# available in source data file: nice_ASVs_long_format.rds
#
nice_asv <- readRDS(file = "./source_data/nice_ASVs_long_format.rds")

## OTUs 16S (MGnify)
nice_otu_raw <- data.table::fread('https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001922/pipelines/5.0/file/ERP024265_taxonomy_abundances_SSU_v5.0.tsv')
#
#select 16S samples (because the file includes 16S and 18S samples)
nice_otu_raw <- nice_otu_raw %>% 
  select(-c(ERR2044671, ERR2044672, ERR2044673,
            ERR2044674, ERR2044675, ERR2044676,
            ERR2044677, ERR2044678, ERR2044679))
#associate colnames with nice samples
#
nice_otu <- 
  nice_otu_raw %>% 
  rename(Taxa = "#SampleID") %>% 
  filter(Taxa != str_detect(Taxa, "sk__Eukaryota")) %>% 
  select(contains("ERR")) %>% 
  mutate(ID = paste0("OTU_", rownames(.))) %>%
  ulrb::prepare_tidy_data(sample_names = colnames(nice_otu_raw[,2:10]))


# mOTU #
nice_mOTU_raw <- data.table::fread('https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001869/pipelines/5.0/file/ERP016733_taxonomy_abundances_SSU_v5.0.tsv')
#
#associate colnames with nice samples
#
nice_mOTU <- 
  nice_mOTU_raw %>% 
  rename(Taxa = "#SampleID") %>% 
  filter(Taxa != str_detect(Taxa, "sk__Eukaryota")) %>% 
  select(contains("ERR")) %>% 
  mutate(ID = paste0("OTU_", rownames(.))) %>%
  ulrb::prepare_tidy_data(sample_names = colnames(nice_mOTU_raw[,2:10]))

## Prepare each dataset for merging
#
nice_asv_clean <- 
  nice_asv %>% 
  filter(!is.na(Kingdom), Family != "Mitochondria", Order != "Chloroplast") %>% 
  mutate(Type = "ASV") %>% 
  select(ID = ASV, Type, Abundance, Sample)

nice_otu <- 
  nice_otu %>% 
  mutate(Type = "OTU") %>%  
  select(Sample, Abundance, Type, ID) 

nice_motu <- 
  nice_mOTU %>% 
  mutate(Type = "mOTU") %>% 
  select(Sample, Abundance, Type, ID) 


# Combined phylogenetic units for N-ICE
nice_dataset <- 
  nice_asv_clean %>% 
  rbind(nice_otu) %>% 
  rbind(nice_motu)


# Standardize sample ID names
nice_dataset <- nice_dataset %>% 
  mutate(Sample = case_when(Sample == "ERR2017139" ~ 'NB_5',
                            Sample == "ERR2017140" ~ 'NB_50',
                            Sample == "ERR2017141" ~ 'NB_250',
                            Sample == "ERR2017142" ~ 'TR_5',
                            Sample == "ERR2017143" ~ 'TR_50',
                            Sample == "ERR2017144" ~ 'TR_250',
                            Sample == "ERR2017145" ~ 'YP_5',
                            Sample == "ERR2017146" ~ 'YP_20',
                            Sample == "ERR2017147" ~ 'YP_250',
                            Sample == "ERR2044662" ~ 'NB_5',
                            Sample == "ERR2044663" ~ 'NB_50',
                            Sample == "ERR2044664" ~ 'NB_250',
                            Sample == "ERR2044665" ~ 'TR_5',
                            Sample == "ERR2044666" ~ 'TR_50',
                            Sample == "ERR2044667" ~ 'TR_250',
                            Sample == "ERR2044668" ~ 'YP_5',
                            Sample == "ERR2044669" ~ 'YP_20',
                            Sample == "ERR2044670" ~ 'YP_250',
                            TRUE ~ Sample))

# Clean before ulrb
## remove singletons
## rarefy for each strategy
asv_rarefaction <- 10000
otu_rarefaction <- asv_rarefaction
motu_rarefaction <- 1000

## note: will give a warning, this is not a problem
set.seed(123)
nice_ulrb_rarefied <- 
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
  mutate(Abundance = Rarefied_abundance, Rarefied_abundance = NULL) %>% 
  group_by(Type) %>% 
  define_rb()

## For comparison with common thresholds approaches
# Make equivalent dataset for other thresholds
nice_ulrb_rarefied_with_thresholds <- 
  nice_ulrb_rarefied %>% 
  group_by(Sample, Type) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(Rarity_0.1_1 = case_when(RelativeAbundance >= 1 ~ "Abundant",
                                  RelativeAbundance <= 0.1 ~ "Rare",
                                  TRUE ~ "Undetermined"),
         Rarity_0.1 = ifelse(RelativeAbundance <= 0.1, "Rare", "Abundant"))

## ulrb vs thresholds
nice_summary_ulrb_rarefied_tidy <- 
  nice_ulrb_rarefied_with_thresholds %>% 
  pivot_longer(cols = c(Classification, Rarity_0.1_1, Rarity_0.1),
               names_to = "Definition", values_to = "Classification") %>% 
  group_by(Sample, Type, Definition, Classification) %>% 
  summarise(Count = n()) %>% 
  mutate(RelativeCount = Count * 100/sum(Count)) %>% 
  mutate(Classification = factor(Classification, 
                                 levels = c("Rare", "Undetermined", "Abundant"))) %>% 
  mutate(Definition = case_when(
    Definition == "Classification" ~ "ulrb",
    Definition == "Rarity_0.1" ~ "one threshold (0.1%)",
    TRUE ~ "two thresholds (0.1% and 1%)")) %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU")))

# Plot for multiple panel figure
gridExtra::grid.arrange(
  #
  nice_ulrb_rarefied %>% 
    group_by(Type, Sample) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
    ungroup() %>% 
    mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
    ggplot(aes(x = reorder(ID, -RelativeAbundance), 
               y = RelativeAbundance, 
               col = Classification)) +
    #geom_point()+
    stat_summary(alpha = 0.5) +
    facet_grid(~Type, scales = "free") + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    scale_y_log10() + 
    geom_hline(yintercept = c(1, 0.1, 0.01), lty = "dashed") +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          axis.ticks.x = element_blank(),
          legend.position = "top") + 
    labs(x = "ranked phylogenetic units",
         y = "mean \U00B1 sd relative abundance (%) \n(Log10 scale)",
         col = "classification: ",
         tag = "a")
  ,
  # silhouette scores
  nice_ulrb_rarefied %>% 
    mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
    ggplot(aes(reorder(ID, -Silhouette_scores), 
               Silhouette_scores, col = Classification)) + 
    stat_summary(alpha = 0.5) + 
    facet_wrap(~Type, scales = "free_x") + 
    #  scale_fill_manual(values = qualitative_colors[c(2,4,7)]) + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          panel.grid = element_blank()) +
    guides(color = "none") + 
    labs(y = "mean \U00B1 sd silhouette score \n(Log10 scale)",
         x = "ranked phylogenetic units",
         tag = "b") + 
    geom_hline(yintercept = 0),
  # alpha diversity
  nice_summary_ulrb_rarefied_tidy %>% 
    ggplot(aes(x = Type, y = Count)) + 
    geom_half_boxplot(side = "l", aes(fill = Classification), 
                      outlier.colour = "red",
                      outlier.shape = "cross",
                      outlier.size = 0.5)+
    geom_half_point(side = "r", size = 0.75, aes(col = Classification)) + 
    facet_grid(~Definition) + 
    scale_fill_manual(values = qualitative_colors[c(3,4,7)]) + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    theme(legend.position = "top",
          strip.background = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey"),
          strip.text = element_text(size = 12)) + 
    labs(x = "kind of phylogenetic unit",
         y = "number of phylogenetic units",
         tag= "c") + 
    guides(color = "none", fill = "none")
)


## report on silhouette scores
# per species
silhouette_strucure_per_observations <- nice_ulrb_rarefied %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
  group_by(Sample, Type, Classification) %>% 
  mutate(Structure = case_when(Silhouette_scores > 0.74 ~ "Strong",
                               Silhouette_scores > 0.50 ~ "Reasonable",
                               Silhouette_scores >= 0 ~ "Weak",
                               Silhouette_scores < 0 ~ "Artificial")) %>% 
  count(Structure) %>%
  ungroup() %>% 
  group_by(Type, Classification, Structure) %>% 
  summarise(total = sum(n)) %>% 
  mutate(ratio = total*100/sum(total)) 

# Plot
plot_sil_species <-
  silhouette_strucure_per_observations %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Classification, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  facet_grid(~Type) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "Percentage of phylogenetic units (%)",
       fill = "Evaluation",
       title = "Position of phylogenetic units in their clusters",
       tag = "a")



# average over each classification
silhouette_structure_per_classifications <- nice_ulrb_rarefied %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
  group_by(Sample, Type, Classification) %>%
  summarise(average_Silhouette_score = mean(Silhouette_scores)) %>% 
  mutate(Structure = case_when(average_Silhouette_score > 0.74 ~ "Strong",
                               average_Silhouette_score > 0.50 ~ "Reasonable",
                               average_Silhouette_score >= 0 ~ "Weak",
                               average_Silhouette_score < 0 ~ "Artificial")) %>% 
  ungroup() %>% 
  group_by(Type, Classification, Structure) %>% 
  count(Structure) %>% 
  summarise(total = sum(n)) %>% 
  mutate(ratio = total*100/sum(total)) 

plot_sil_class <- silhouette_structure_per_classifications %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Classification, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  facet_grid(~Type) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "Percentage of samples (%)",
       fill = "Evaluation",
       title = "Structure of each classification",
       tag = "b")

# average overall (i.e. structure of the clustering)
silhouette_structure_overall <- nice_ulrb_rarefied %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
  group_by(Sample, Type) %>% 
  summarise(average_Silhouette_score = mean(Silhouette_scores)) %>% 
  mutate(Structure = case_when(average_Silhouette_score > 0.74 ~ "Strong",
                               average_Silhouette_score > 0.50 ~ "Reasonable",
                               average_Silhouette_score >= 0 ~ "Weak",
                               average_Silhouette_score < 0 ~ "Artificial")) %>% 
  ungroup() %>% 
  group_by(Type, Structure) %>% 
  count(Structure) %>% 
  summarise(total = sum(n)) %>% 
  mutate(ratio = total*100/sum(total)) 

#
plot_sil_nice_overall <- silhouette_structure_overall %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Type, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "Percentage of samples (%)",
       fill = "Evaluation",
       title = "Structure of final cluster",
       tag = "c")


### overall silhouette evaluation
grid.arrange(
  plot_sil_species + theme(axis.text.x = element_text(angle = 90, vjust  =0)),
  arrangeGrob(
    plot_sil_class + theme(axis.text.x = element_text(angle = 90, vjust  =0)) + guides(fill = "none"),
    arrangeGrob(
      plot_sil_nice_overall + guides(fill = "none"),
      ncol = 1)),
  nrow = 1, ncol = 2)


## why are ASVs harder to cluster?

nice_ulrb_rarefied %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>%
  ggplot(aes(Abundance, fill = Type)) + 
  geom_histogram(bins = 30, alpha = 0.75) +
  facet_wrap(~Type, scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values = qualitative_colors) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top")

##

## phylogenetic units benchmark
set.seed(123)
nice_ulrb_rarefied_for_bench <- 
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
  mutate(Abundance = Rarefied_abundance, Rarefied_abundance = NULL)

#
nice_ulrb_rarefied_for_bench_ASV <- 
  nice_ulrb_rarefied_for_bench %>% filter(Type == "ASV")
nice_ulrb_rarefied_for_bench_OTU <- 
  nice_ulrb_rarefied_for_bench %>% filter(Type == "OTU")
nice_ulrb_rarefied_for_bench_mOTU <- 
  nice_ulrb_rarefied_for_bench %>% filter(Type == "mOTU")

#
phylgeneticUnits_benchmark <- 
  microbenchmark(ASV = {define_rb(nice_ulrb_rarefied_for_bench_ASV)},
                 OTU = {define_rb(nice_ulrb_rarefied_for_bench_OTU)},
                 mOTU = {define_rb(nice_ulrb_rarefied_for_bench_mOTU)})
#
phylgeneticUnits_benchmark %>% 
  autoplot() + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) + 
  labs(x = "time (miliseconds)",
       y = "Phylogenetic unit",
       title = "Time performance of define_rb()",
       subtitle = "100 replications")


## Sequencing summary statistics

