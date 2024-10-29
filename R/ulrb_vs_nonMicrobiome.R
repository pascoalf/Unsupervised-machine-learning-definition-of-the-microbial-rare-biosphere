# ulrb vs non-microbiome datasets

## BCI datasets, see "BCI_processing.R" script.
#BCI_tidy <- prepare_tidy_data(BCI, samples_in = "rows", sample_names = rownames(BCI))
# 
#BCI_ulrb <- BCI_tidy %>%  define_rb()

## Ants dataset
data("antsA")
ants_sample_names <- rownames(antsA)
antsA_tidy <- antsA %>% 
  prepare_tidy_data(samples_in = "rows", sample_names = ants_sample_names) %>% 
  filter(Abundance > 0) %>% 
  filter(Sample != "sam_95")
#
antsA_tidy %>%
  group_by(Sample) %>% 
  distinct(Abundance) %>% 
  count() %>% 
  pull(n) %>% 
  min()
#
ants_ulrb <- define_rb(antsA_tidy)
# rank ant species
ranked_ants <- antsA_tidy %>% 
  group_by(Taxa_id) %>% 
  summarise(totalAbundance = sum(Abundance)) %>% 
  arrange(desc(totalAbundance)) %>% 
  pull(Taxa_id)


# Ants RAC
Ants_rac <- ants_ulrb %>% 
    group_by(Sample) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>%
    ungroup() %>% 
    mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
    group_by(Sample) %>% 
    arrange(desc(Abundance)) %>% 
    mutate(uniqueRank = row_number()) %>% 
    ungroup() %>% 
    mutate(Taxa_id = factor(Taxa_id, levels = ranked_ants)) %>% 
    ggplot(aes(uniqueRank, y = RelativeAbundance, col = Classification)) +
    geom_line(aes(group = Group), col = "grey72") + 
    geom_jitter(alpha = 0.75, size = 3, width = 0.1) +
    #stat_summary() +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(), 
          strip.background = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12)) + 
    labs(x = "Ranked ant species",
         y = "Relative abundance (%)",
         col = "Classification:",
         title = "Ants dataset",
         tag = "a") +
    scale_color_manual(values = qualitative_colors[c(3,5,7)])

# Ants silhouette
Ants_sil <- ants_ulrb %>% 
  mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
  group_by(Group) %>% 
  arrange(desc(Silhouette_scores)) %>% 
  mutate(uniqueRank = row_number()) %>% 
  ungroup() %>%
  ggplot(aes(uniqueRank, 
             Silhouette_scores,
             col = Classification)) +
  geom_hline(yintercept = 0)  + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  geom_line(aes(group = Group), col = "grey72")+
  geom_jitter(alpha = 0.75,size = 3, width = 0.1) + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  labs(x = "Ranked ant species",
       y = "Average Silhouette score",
       col = "Classification:",
       title = "Ants dataset",
       tag = "b") +
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

# BCI - ulrb
BCI_RAC <- BCI_ulrb %>%  
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
  group_by(Group) %>% 
  arrange(desc(Abundance)) %>% 
  mutate(uniqueRank = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(uniqueRank,
             RelativeAbundance,
             col = Classification)) + 
  geom_line(aes(group = Group), col = "grey72") + 
  geom_point(alpha = 0.75, size = 1.5) +
  theme_classic() + 
  scale_y_log10() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  labs(x = "Ranked ant species",
       y = "Relative abundance (%)",
       col = "Classification:",
       title = "BCI dataset",
       tag = "c") +
  scale_color_manual(values = qualitative_colors[c(3, 5, 7)])

# BCI silhouette
BCI_sil <- BCI_ulrb %>% 
  mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
  group_by(Group) %>% 
  arrange(desc(Silhouette_scores)) %>% 
  mutate(uniqueRank = row_number()) %>% 
  ungroup() %>%
  ggplot(aes(uniqueRank,
             Silhouette_scores,
             col = Classification)) +
  geom_line(aes(group = Group), col = "grey72") + 
  geom_point(size = 1.5, alpha = 0.75) + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  labs(title = "BCI dataset",
       x = "Ranked tree species",
       y = "Average Silhouette score",
       col = "Classification:",
       tag = "d") +
  scale_color_manual(values = qualitative_colors[c(3, 5, 7)]) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(0.5), lty = "dashed")

#  Plot BCI and Plants at once
grid.arrange(
Ants_rac,
Ants_sil,
BCI_RAC,
BCI_sil)


### more tests
ants_absolute <- ants_ulrb
ants_relative <- antsA_tidy %>% 
  group_by(Sample) %>% 
  mutate(Abundance = Abundance*100/sum(Abundance)) %>% 
  define_rb()

classifications_absolute_values <- ants_absolute %>% select(Sample, Classification, Taxa_id) %>% mutate(Type = "Absolute")
classifications_relative_values <- ants_relative %>% select(Sample, Classification, Taxa_id) %>% mutate(Type = "Relative")

comparison <- classifications_absolute_values == classifications_relative_values

combined_ants <- classifications_absolute_values %>% bind_rows(classifications_relative_values)


BCI_absolute <- BCI_ulrb
BCI_relative <- BCI_tidy %>% 
  group_by(Sample) %>% 
  mutate(Abundance = Abundance*100/sum(Abundance)) %>% 
  define_rb()

#
classifications_bci_abs <- BCI_absolute %>% select(Sample, Classification, Taxa_id) %>% mutate(Type = "Absolute")
classifications_bci_relative <- BCI_relative %>% select(Sample, Classification, Taxa_id) %>% mutate(Type = "Relative")

combined_classifications_bci <- classifications_bci_abs %>% rbind(classifications_bci_relative)


gridExtra::grid.arrange(
combined_ants %>% 
  group_by(Sample, Type) %>% 
  count(Classification)  %>% 
  ggplot(aes(Classification, n, col = Classification)) +
  geom_point() + 
  facet_grid(~Type) + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  scale_color_manual(values = qualitative_colors[c(3, 5, 7)]) + 
  labs(title = "Ants dataset",
       x = "Classification",
       y = "Count",
       col = "Classification:",
       tag = "a") ,
#
combined_classifications_bci %>% 
  group_by(Sample, Type) %>% 
  count(Classification)  %>% 
  ggplot(aes(Classification, n, col = Classification)) +
  geom_point() + 
  facet_grid(~Type) + 
  theme_classic() + 
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  scale_color_manual(values = qualitative_colors[c(3, 5, 7)]) + 
  labs(title = "BCI dataset",
       x = "Classification",
       y = "Count",
       col = "Classification:",
       tag = "b"))


