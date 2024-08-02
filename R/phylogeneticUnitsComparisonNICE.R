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


## Add FuzzyQ classification
# ASVs
nice_ASVs_matrix <- nice_asv_clean %>% 
  ungroup() %>%   
  select(Sample, Abundance, ID, Type) %>% 
  mutate(Sample_unique = paste(Sample, Type)) %>%
  select(Sample_unique, Abundance, ID) %>% 
  tidyr::pivot_wider(names_from = "ID",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())
#
rownames(nice_ASVs_matrix) <- nice_ASVs_matrix$Sample_unique
nice_ASVs_matrix$Sample_unique <- NULL
nice_ASVs_matrix <- as.matrix(nice_ASVs_matrix)

#
ASVs_fuzzyQ <- fuzzyq(nice_ASVs_matrix, sorting = TRUE)$spp
ASVs_fuzzyQ$Type = "ASV"

# OTUs
nice_otu_matrix <- nice_otu %>% 
  ungroup() %>%   
  select(Sample, Abundance, ID) %>% 
  tidyr::pivot_wider(names_from = "ID",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())
#
rownames(nice_otu_matrix) <- nice_otu_matrix$Sample
nice_otu_matrix$Sample <- NULL
nice_otu_matrix <- as.matrix(nice_otu_matrix)

#
OTUs_fuzzyQ <- fuzzyq(nice_otu_matrix, sorting = TRUE)$spp
OTUs_fuzzyQ$Type = "OTU"

# mOTUs
nice_mOTU_matrix <- nice_motu %>% 
  ungroup() %>%   
  select(Sample, Abundance, ID) %>% 
  tidyr::pivot_wider(names_from = "ID",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())
#
rownames(nice_mOTU_matrix) <- nice_mOTU_matrix$Sample
nice_mOTU_matrix$Sample <- NULL
nice_mOTU_matrix <- as.matrix(nice_mOTU_matrix)

#
mOTUs_fuzzyQ <- fuzzyq(nice_mOTU_matrix, sorting = TRUE)$spp
mOTUs_fuzzyQ$Type = "mOTU"

## combine fuzzyQ results
nice_fuzyQ <- rbind(ASVs_fuzzyQ, 
      OTUs_fuzzyQ, 
      mOTUs_fuzzyQ) %>% 
  mutate(Classification = ifelse(cluster == 0, "Rare", "Common"),
         Definition = "FuzzyQ",
         ID = rownames(.))

#
nice_fuzyQ_full <- nice_dataset %>% 
  left_join(nice_fuzyQ)

#
nice_fuzzy_diversity <- nice_fuzyQ_full %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample, Type, Classification) %>% 
  summarise(Count = specnumber(Abundance), 
            Definition = "FuzzyQ") 

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
    geom_point(alpha = 0.5, size = 2)+
    #stat_summary(alpha = 0.5) +
    facet_grid(~Type, scales = "free") + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    scale_y_log10() + 
    geom_hline(yintercept = c(1, 0.1, 0.01), lty = "dashed") +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12)) + 
    labs(x = "Ranked phylogenetic units",
         y = "Relative abundance (%) \n(Log10 scale)",
         col = "Classification: ",
         tag = "a")
  ,
  # silhouette scores
  nice_ulrb_rarefied %>% 
    mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
    ggplot(aes(reorder(ID, -Silhouette_scores), 
               Silhouette_scores, col = Classification)) + 
    geom_point(alpha = 0.5, size = 2) + 
    #stat_summary(alpha = 0.5) + 
    facet_wrap(~Type, scales = "free_x") + 
    #  scale_fill_manual(values = qualitative_colors[c(2,4,7)]) + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14),
          panel.grid = element_blank(),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12)) +
    guides(color = "none") + 
    labs(y = "Silhouette score \n(Log10 scale)",
         x = "Ranked phylogenetic units",
         tag = "b") + 
    geom_hline(yintercept = c(0)) + 
    geom_hline(yintercept = c(-0.5, 0.5), lty = "dashed", color = "grey41")
)

##
# alpha diversity
nice_summary_ulrb_rarefied_tidy %>% 
  filter(!is.na(Classification)) %>% 
  mutate(Classification = factor(Classification, 
                                 levels = c("Rare", "Undetermined", "Abundant", "Common"))) %>%
  mutate(Definition = factor(Definition, 
                             levels = c("one threshold (0.1%)",
                                        "two thresholds (0.1% and 1%)",
                                        "ulrb"))) %>% 
  ggplot(aes(x = Classification, y = Count)) +
  geom_half_boxplot(side = "l", aes(fill = Type), 
                    outlier.colour = "red",
                    outlier.shape = "cross",
                    outlier.size = 2)+
  geom_half_point(side = "r", size = 1.5, aes(col = Type)) +
  stat_summary(aes(y = Count, group = Type, 
                   color = Type), 
               fun = median, geom = "line",
               lwd = 0.5) +
  facet_grid(~Definition, scales = "free_x") + 
  scale_fill_manual(values = qualitative_colors[c(1, 2, 3)]) + 
  scale_color_manual(values = qualitative_colors[c(1, 2, 3)]) + 
  theme_bw() + 
  theme(legend.position = "top",
        strip.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)) + 
  labs(x = "Classification",
       y = "Number of phylogenetic units",
       fill = "Phylogenetic units: ",
       col = "Phylogenetic units: ")


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
       fill = "Evaluation: ",
       title = "Evaluation of unit within its classification/cluster",
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
       title = "Evaluation of each classification/cluster",
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
       title = "Evaluation of final cluster \n(by average Silhouette score)",
       tag = "c")


### overall silhouette evaluation
grid.arrange(
  plot_sil_species + theme(axis.text.x = element_text(angle = 90, vjust = 0)),
  arrangeGrob(
    plot_sil_class + theme(axis.text.x = element_text(angle = 90, vjust = 0)) + 
      guides(fill = "none"),
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
  scale_x_log10() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  labs(y = "Count",
       x = "Abundance (Log10 scale)",
       fill = "Phylogenetic unit:")

## different values is a simpler way 
nice_ulrb_rarefied %>% 
  group_by(Type) %>% 
  count(Abundance) %>% 
  summarise(n(), min(Abundance), max(Abundance))

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
#phylgeneticUnits_benchmark <- 
#  microbenchmark(ASV = {define_rb(nice_ulrb_rarefied_for_bench_ASV)},
#                 OTU = {define_rb(nice_ulrb_rarefied_for_bench_OTU)},
#                 mOTU = {define_rb(nice_ulrb_rarefied_for_bench_mOTU)})
#
#phylgeneticUnits_benchmark %>% 
#  autoplot() + 
#  theme_bw() +
#  theme(panel.grid.minor = element_blank(),
#        panel.grid.major.y = element_blank()) + 
#  labs(x = "time (miliseconds)",
#       y = "Phylogenetic unit",
#       title = "Time it takes to run define_rb()",
#       subtitle = "100 replications")



#####
gridExtra::grid.arrange(
nice_fuzyQ %>% 
  ggplot(aes(reorder(ID, -Common.I), 
             y = Common.I, 
             col = Classification)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  labs(y = "Commonnality index",
       x = "Ranked pylogenetic units",
       fill = "Classification: ") + 
  facet_grid(~Type, scales = "free_x") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  ylim(c(0,1)),

# Silhouette scores
nice_fuzyQ %>% 
  ggplot(aes(reorder(ID, -sil_width),
             y = sil_width, 
             fill = Classification,
             col = Classification)) + 
  geom_col() + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size =14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  facet_grid(~Type, scales = "free_x") + 
  labs(y = "Silhouette score",
       x = "Phylogenetic unit")+ 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]) + 
  guides(fill = FALSE, col = FALSE)
)
  



