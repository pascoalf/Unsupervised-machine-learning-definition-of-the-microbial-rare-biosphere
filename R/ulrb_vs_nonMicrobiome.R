# ulrb vs non-microbiome datasets

## BCI datasets
# transform to long format
BCI_tidy <- prepare_tidy_data(BCI, samples_in = "rows", sample_names = rownames(BCI))
# 
BCI_ulrb <- BCI_tidy %>% 
  define_rb()

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
