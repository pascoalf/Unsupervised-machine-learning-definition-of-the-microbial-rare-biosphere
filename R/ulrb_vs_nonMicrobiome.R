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



#  Plot BCI and Plants at once

grid.arrange(
  ants_ulrb %>%
    group_by(Sample) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>%
    ungroup() %>% 
    mutate(Taxa_id = factor(Taxa_id, levels = ranked_ants)) %>% 
    ggplot(aes(Taxa_id, y = Abundance, col = Classification)) +
    #geom_point(size = 1.5) +
    stat_summary() +
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
         y = "Mean \U00B1 sd\n relative abundance (%)",
         col = "Classification:",
         title = "Ants dataset",
         tag = "a") +
    scale_color_manual(values = qualitative_colors[c(3,4,7)]),
  ants_ulrb %>% 
    ggplot(aes(reorder(Taxa_id, -Silhouette_scores), 
               Silhouette_scores,
               col = Classification)) + 
    stat_summary() +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(), 
          strip.background = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12)) + 
    geom_hline(yintercept = 0)  + 
    geom_hline(yintercept = 0.5, lty = "dashed") + 
    labs(x = "Ranked ant species",
         y = "Mean \U00B1 sd of\n average Silhouette score",
         col = "Classification:",
         title = "Ants dataset",
         tag = "b") +
    scale_color_manual(values = qualitative_colors[c(3,4,7)]),
# BCI - ulrb
  BCI_ulrb %>% 
    group_by(Sample) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
    ungroup() %>% 
    ggplot(aes(reorder(Taxa_id, -RelativeAbundance),
               RelativeAbundance,
               col = Classification)) + 
    stat_summary() + 
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
         y = "Mean \U00B1 sd \nrelative abundance (%)",
         col = "Classification:",
         title = "BCI dataset",
         tag = "c") +
    scale_color_manual(values = qualitative_colors[c(3, 4, 7)]),
  BCI_ulrb %>% 
    ggplot(aes(reorder(Taxa_id, -Silhouette_scores),
               Silhouette_scores,
               col = Classification)) + 
    stat_summary() + 
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
         y = "Mean \U00B1 sd\n of average Silhouette score",
         col = "Classification:",
         tag = "d") +
    scale_color_manual(values = qualitative_colors[c(3, 4, 7)]) + 
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = c(0.5), lty = "dashed"))
