## MOSJ 2019
# Short vs full 16S rRNA gene (section 2)

## Clean and merge taxonomy with abundance tables
# Short 16S
ASV_abundance_short <- readRDS("./source_data/seqtab.nochim_mosj.rds") %>% 
  t() %>% 
  as.data.frame() 

ASV_abundance_short$Sequence <- rownames(ASV_abundance_short)
rownames(ASV_abundance_short) <- NULL
#

# Silva taxonomy
ASV_taxonomy_silva <- readRDS("./source_data/short_mosj_taxa_silva.rds") %>% as.data.frame()
ASV_taxonomy_silva$Sequence <- rownames(ASV_taxonomy_silva)
rownames(ASV_taxonomy_silva) <- NULL

#
ASV_taxonomy_silva <- readRDS("./source_data/short_mosj_taxa_silva.rds") %>% as.data.frame() %>% 
  mutate(Sequencing_technology = "Illumina",
         Marker = "V4-V5 16S rRNA gene",
         Database = "Silva",
         Sequence = rownames(.)) %>% 
  select(-Genus)

rownames(ASV_taxonomy_silva) <- NULL

short_asv <- ASV_taxonomy_silva %>% 
  left_join(ASV_abundance_short, by = "Sequence")

# Full-length 16S
ASV_abundance_full_mosj <- 
  readRDS("./source_data/seq_tab_nochim_full_mosj.rds") %>% 
  t() %>% as.data.frame()
ASV_abundance_full_mosj$Sequence <- rownames(ASV_abundance_full_mosj)
rownames(ASV_abundance_full_mosj) <- NULL

# Silva taxonomy (full)
ASV_taxonomy_silva_full_mosj_df <- readRDS("./source_data/ASV_taxonomy_silva_full_mosj.rds") %>% as.data.frame()
ASV_taxonomy_silva_full_mosj_df$Sequence <- rownames(ASV_taxonomy_silva_full_mosj_df)
rownames(ASV_taxonomy_silva_full_mosj_df) <- NULL

ASV_abundance_taxonomy_full_mosj_df_silva <- ASV_abundance_full_mosj %>% 
  left_join(ASV_taxonomy_silva_full_mosj_df, by = "Sequence") %>% 
  mutate(Sequencing_technology = "PacBio",
         Marker = "Full-length 16S rRNA gene",
         Database = "Silva") %>% 
  mutate(Species = ifelse(is.na(Species), NA, paste(Genus, Species)))

# Select relevant cols and clean taxonomy if necessary
short_asv_clean <- short_asv %>% 
  filter(!is.na(Kingdom), 
         Kingdom != "Eukaryota",
         Order != "Chlorpolast",
         Family != "Mitochondria") %>% 
  mutate(ID = paste0("ASV_", row_number(.))) %>% 
  pivot_longer(cols = contains("M19"), values_to = "Abundance", names_to = "Sample")
#
short_asv_df <- short_asv_clean %>% 
  select(ID, Sample, Marker, Abundance) %>% 
  mutate(Sample = str_remove(Sample, "_R1_001_filt.fastq.gz"),
         Sample = str_replace(Sample, "-", "_"))

# clean full length table
full_asv_df <- 
  ASV_abundance_taxonomy_full_mosj_df_silva %>% 
  mutate(ID = paste0("ASV_", row_number())) %>% 
  filter(Order != "Chloroplast", Family != "Mitochondria") %>% 
  pivot_longer(cols = contains("M19"), values_to = "Abundance", names_to = "Sample") %>%
  select(ID, Sample, Marker, Abundance) %>% 
  mutate(Sample = str_remove(Sample, ".hifi_reads.fastq"))


## Combine short and full
short_and_full_asv <- rbind(short_asv_df, full_asv_df)

#### Analysis #####
short_asv_rarefaction <- 10000
full_asv_rarefaction <- 5000

# Remove samples that don't satisfy minimum rarefaction threshold of 10 000 reads
# in both strategies at same time.
short_vs_full_samples_to_remove <- 
  short_and_full_asv %>% 
  group_by(Marker, Sample) %>% 
  summarise(total_reads = sum(Abundance)) %>% 
  filter(total_reads < 10000) %>% 
  pull(Sample)

# Define rare biosphere
set.seed(123)
# will give a warning, this is not a problem
short_vs_full_ulrb <- 
  short_and_full_asv %>% 
  filter(!Sample %in% short_vs_full_samples_to_remove) %>% 
  mutate(Rarefaction = case_when(str_detect(Marker, "V4-V5") ~ short_asv_rarefaction,
                                 str_detect(Marker, "Full") ~ full_asv_rarefaction)) %>% 
  group_by(Marker, Sample, Rarefaction) %>% 
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(.x$Abundance, 
                                          sample = Rarefaction))))) %>% 
  unnest(c(data, Rarefied_reads)) %>% 
  rename(Rarefied_abundance = "V1") %>% 
  mutate(Abundance = Rarefied_abundance, Rarefied_abundance = NULL) %>% 
  group_by(Marker, Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>% ## abundance to relative abundance
  ungroup() %>% 
  group_by(Marker) %>% 
  define_rb()

## Compare with other definitions
short_vs_full_ulrb_with_thresholds <- 
  short_vs_full_ulrb %>% 
  group_by(Sample, Marker) %>% 
  ungroup() %>% 
  mutate(Rarity_0.1_1 = case_when(RelativeAbundance >= 1 ~ "Abundant",
                                  RelativeAbundance <= 0.1 ~ "Rare",
                                  TRUE ~ "Undetermined"),
         Rarity_0.1 = ifelse(RelativeAbundance <= 0.1, "Rare", "Abundant"))

short_vs_full_ulrb_with_thresholds_summary <- 
  short_vs_full_ulrb_with_thresholds %>% 
  pivot_longer(cols = c(Classification, Rarity_0.1_1, Rarity_0.1),
               names_to = "Definition", values_to = "Classification") %>% 
  group_by(Sample, Marker, Definition, Classification) %>% 
  summarise(Count = n()) %>% 
  mutate(RelativeCount = Count * 100/sum(Count)) %>% 
  mutate(Classification = factor(Classification, 
                                 levels = c("Rare", "Undetermined", "Abundant"))) %>% 
  mutate(Definition = case_when(Definition == "Classification" ~ "ulrb",
                                Definition == "Rarity_0.1" ~ "one threshold (0.1%)",
                                TRUE ~ "two thresholds (0.1% and 1%)"))
# View results
gridExtra::grid.arrange(
  short_vs_full_ulrb %>%
    ggplot(aes(reorder(ID, -RelativeAbundance), 
               RelativeAbundance, 
               col = Classification)) + 
    stat_summary(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = c(0.01, 0.1, 1), linetype = "dashed") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x =element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          strip.text = element_text(size = 10)) + 
    facet_grid(~Marker, scales = "free") + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    scale_y_log10() + 
    labs(x = "ranked ASV",
         y = "mean \U00B1 sd relative abundance (%) \n(Log10 scale)",
         col = "classification: ")
  ,
  # Silhouette plot
  short_vs_full_ulrb %>%
    ggplot(aes(reorder(ID, -Silhouette_scores), 
               Silhouette_scores, 
               col = Classification)) + 
    stat_summary(size = 0.5, alpha = 0.5) +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x =element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          strip.text = element_text(size = 10)) + 
    facet_grid(~Marker, scales = "free") + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    guides(color = "none") + 
    labs(x = "ranked ASV",
         y = "mean \U00B1 sd silhouette score") +
    geom_hline(yintercept = 0), 
  short_vs_full_ulrb_with_thresholds_summary %>% 
    ggplot(aes(Marker, Count)) + 
    facet_grid(~Definition) + 
    geom_half_boxplot(side = "l", aes(fill = Classification),
                      outlier.colour = "red",
                      outlier.shape = "cross",
                      outlier.size = 0.5) +
    geom_half_point(side = "r", aes(col = Classification)) + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    scale_fill_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() + 
    theme(legend.position = "top",
          strip.background = element_blank(),
          panel.grid.major.x = element_line(color = "grey"),
          panel.grid.minor.y = element_blank()) + 
    labs(x = "marker gene",
         y = "number of ASVs") + 
    guides(fill = "none", color = "none")
)

### deeper view on silhouette scores
# per species
(short_vs_full_silhouette_species <- 
    short_vs_full_ulrb %>% 
    group_by(Sample, Marker, Classification) %>% 
    mutate(Structure = case_when(Silhouette_scores > 0.74 ~ "Strong",
                                 Silhouette_scores > 0.50 ~ "Reasonable",
                                 Silhouette_scores >= 0 ~ "Weak",
                                 Silhouette_scores < 0 ~ "Artificial")) %>% 
    count(Structure) %>%
    ungroup() %>% 
    group_by(Marker, Classification, Structure) %>% 
    summarise(total = sum(n)) %>% 
    mutate(ratio = total*100/sum(total))) 

#
short_vs_full_plot_sil_species <-
  short_vs_full_silhouette_species %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Classification, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  facet_grid(~Marker) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "Percentage of phylogenetic units (%)",
       fill = "Evaluation",
       title = "Positioning of phylogenetic units in their clusters",
       tag = "a")

##
# average over each classification
short_vs_full_silhouette_classifications <- 
  short_vs_full_ulrb %>% 
  group_by(Sample, Marker, Classification) %>%
  summarise(average_Silhouette_score = mean(Silhouette_scores)) %>% 
  mutate(Structure = case_when(average_Silhouette_score > 0.74 ~ "Strong",
                               average_Silhouette_score > 0.50 ~ "Reasonable",
                               average_Silhouette_score >= 0 ~ "Weak",
                               average_Silhouette_score < 0 ~ "Artificial")) %>% 
  ungroup() %>% 
  group_by(Marker, Classification, Structure) %>% 
  count(Structure) %>% 
  summarise(total = sum(n)) %>% 
  mutate(ratio = total*100/sum(total)) 

##
short_vs_full_plot_sil_class <- 
  short_vs_full_silhouette_classifications %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Classification, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  facet_grid(~Marker) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "Percentage of samples (%)",
       fill = "Evaluation",
       title = "Structure of each classification",
       tag = "b")

##
# average overall (i.e. structure of the clustering)
short_vs_full_silhouette_structure_overall <- 
  short_vs_full_ulrb %>% 
  group_by(Sample, Marker) %>% 
  summarise(average_Silhouette_score = mean(Silhouette_scores)) %>% 
  mutate(Structure = case_when(average_Silhouette_score > 0.74 ~ "Strong",
                               average_Silhouette_score > 0.50 ~ "Reasonable",
                               average_Silhouette_score >= 0 ~ "Weak",
                               average_Silhouette_score < 0 ~ "Artificial")) %>% 
  ungroup() %>% 
  group_by(Marker, Structure) %>% 
  count(Structure) %>% 
  summarise(total = sum(n)) %>% 
  mutate(ratio = total*100/sum(total)) 

#
plot_short_vs_full_silhouette_structure_overall <- 
  short_vs_full_silhouette_structure_overall %>% 
  mutate(Structure = factor(Structure, c("Strong",
                                         "Reasonable",
                                         "Weak",
                                         "Artificial"))) %>% 
  ggplot(aes(x = Marker, y = ratio, fill = Structure)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "top") + 
  labs(y = "percentage of samples (%)",
       fill = "Evaluation",
       title = "Structure of final clustering",
       tag = "c")

##

grid.arrange(
  short_vs_full_plot_sil_species + theme(axis.text.x = element_text(angle = 90, vjust  = 0)),
  arrangeGrob(
    short_vs_full_plot_sil_class + theme(axis.text.x = element_text(angle = 90, vjust  = 0)) + guides(fill = "none"),
    plot_short_vs_full_silhouette_structure_overall + guides(fill = "none")),
  nrow = 1, ncol = 2)


# Abundance distribution
grid.arrange(
  short_vs_full_ulrb %>% 
    ggplot(aes(Abundance, fill = Marker)) + 
    geom_histogram(bins = 50, alpha = 0.75) + 
    facet_wrap(~Marker, scales = "free_x") +
    theme_bw() + 
    scale_fill_manual(values = qualitative_colors) + 
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10),
          panel.grid = element_blank(),
          legend.position = "top") + 
    labs(tag = "a")
  ,
  short_vs_full_ulrb %>%
    ggplot(aes(Abundance, fill = Classification)) + 
    geom_histogram(bins = 50, alpha = 0.75) + 
    facet_wrap(~Marker, scales = "free_x") +
    theme_bw() + 
    scale_fill_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 10),
          panel.grid = element_blank(),
          legend.position = "top") + 
    labs(tag = "b")
)

# short vs full microbenchmarking #
short_vs_full_for_benchmark <- 
  short_and_full_asv %>% 
  filter(!Sample %in% short_vs_full_samples_to_remove) %>% 
  mutate(Rarefaction = case_when(str_detect(Marker, "V4-V5") ~ short_asv_rarefaction,
                                 str_detect(Marker, "Full") ~ full_asv_rarefaction)) %>% 
  group_by(Marker, Sample, Rarefaction) %>% 
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(.x$Abundance, 
                                          sample = Rarefaction))))) %>% 
  unnest(c(data, Rarefied_reads)) %>% 
  rename(Rarefied_abundance = "V1") %>% 
  mutate(Abundance = Rarefied_abundance, Rarefied_abundance = NULL) %>% 
  group_by(Marker, Sample) %>% 
  mutate(RelativeAbundance = Abundance * 100/sum(Abundance)) %>% ## abundance to relative abundance
  ungroup()

#
short_vs_full_for_benchmark_V4V5 <- 
  short_vs_full_for_benchmark %>% 
  filter(str_detect(Marker, "V4-V5"))
#
short_vs_full_for_benchmark_Full <- 
  short_vs_full_for_benchmark %>% 
  filter(str_detect(Marker, "Full"))

# 
short_vs_full_microbenchmark <- 
  microbenchmark("V4V5 16S rRNA gene" = {define_rb(short_vs_full_for_benchmark_V4V5)},
                 "Full-length 16S rRNA gene" = {define_rb(short_vs_full_for_benchmark_Full)})


short_vs_full_microbenchmark %>% 
  autoplot() + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) + 
  labs(x = "time (miliseconds)",
       y = "Marker gene",
       title = "Time performance of define_rb()",
       subtitle = "100 replications")

