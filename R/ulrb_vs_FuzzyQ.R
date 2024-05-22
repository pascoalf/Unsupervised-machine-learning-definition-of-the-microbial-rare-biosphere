## FuzzyQ vs ulrb ##
data(antsA)

FQAnts <- fuzzyq(antsA[-95,], sorting = TRUE)

# ants fuzzyQ figure
FQAnts_df <- FQAnts$spp %>% 
  mutate(Species = rownames(.),
         Classification = ifelse(cluster == 0, "Rare", "Common"))

# Figure
gridExtra::grid.arrange(
FQAnts_df %>% 
  ggplot(aes(reorder(Species, -Common.I),
             Common.I,
             col = Classification)) + 
  geom_point() + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Commonality index",
       x = "Species",
       title = "FuzzyQ - Ants dataset") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]) + 
  ylim(c(0,1)),
FQAnts_df %>% 
  ggplot(aes(reorder(Species, -sil_width),
             sil_width,
             fill = Classification,
             col = Classification)) + 
  geom_col() + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Silhouette score",
       x = "Species") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]))

# evaluate fuzzyQ - Ants
fuzzy_ants_quality <- 
  FQAnts_df %>% 
  mutate(evaluate_Sil = case_when(sil_width > 0.71 ~ "Strong cluster",
                                  sil_width > 0.51 ~ "Reasonable cluster",
                                  sil_width > 0.26 ~ "Weak cluster",
                                  sil_width <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "FuzzyQ quality - Ants dataset")
  

#
all_years_ulrb_matrix <- all_years_ulrb %>% 
  ungroup() %>%   
  select(Sample, Abundance, Sequence, Depth, year) %>% 
  mutate(Sample_unique = paste(Sample, year)) %>%
  select(Sample_unique, Abundance, Sequence) %>% 
  tidyr::pivot_wider(names_from = "Sequence",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())


rownames(all_years_ulrb_matrix) <- all_years_ulrb_matrix$Sample_unique
all_years_ulrb_matrix$Sample_unique <- NULL
all_years_ulrb_matrix <- as.matrix(all_years_ulrb_matrix)

all_years_fuzzy <- fuzzyq(all_years_ulrb_matrix, sorting = TRUE)


all_years_fuzzy_df <- all_years_fuzzy$spp %>% 
  mutate(ID = rownames(.),
         Classification = ifelse(cluster == 0, "Rare", "Common"))
#
gridExtra::grid.arrange(
all_years_fuzzy_df %>% 
  ggplot(aes(reorder(ID, -Common.I),
             Common.I,
             col = Classification,
             fill = Classification)) + 
  geom_point() + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Commonality index",
       x = "ASVs",
       title = "FuzzyQ - MOSJ2016-2020 dataset") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]) + 
  ylim(c(0,1)),
all_years_fuzzy_df %>% 
  ggplot(aes(reorder(ID, -sil_width),
             sil_width,
             col = Classification,
             fill = Classification)) + 
  geom_col() + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Silhouette score",
       x = "ASVs") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]))

# FuzzyQ quality - mosj2016-2020 dataset
fuzzyQ_quality_mosj <- all_years_fuzzy_df %>%   
  mutate(evaluate_Sil = case_when(sil_width > 0.71 ~ "Strong cluster",
                                  sil_width > 0.51 ~ "Reasonable cluster",
                                  sil_width > 0.26 ~ "Weak cluster",
                                  sil_width <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "FuzzyQ quality - MOSJ2016-2020 dataset")

## fuzzyQ - BCI dataset
# load BCI dataset
data(BCI)
#
FuzzyQ_plants <- fuzzyq(BCI, sorting = TRUE)
FuzzyQ_plants <-  FuzzyQ_plants$spp %>% 
  as.data.frame() %>% 
  mutate(Species = rownames(.),
         Classification = ifelse(cluster == 0, "Rare", "Common"))
#
gridExtra::grid.arrange(
FuzzyQ_plants %>% 
  ggplot(aes(reorder(Species, -Common.I),
             Common.I, 
             col = Classification)) + 
  geom_point() +  
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Commonality index",
       x = "Species",
       title = "FuzzyQ - BCI dataset") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]) + 
  ylim(c(0,1)),
FuzzyQ_plants %>% 
  ggplot(aes(reorder(Species, -sil_width),
             sil_width, 
             col = Classification,
             fill = Classification)) + 
  geom_col() +  
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  theme_classic() + 
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) + 
  labs(y = "Silhouette score",
       x = "Species") + 
  scale_color_manual(values = qualitative_colors[c(1,2)]) + 
  scale_fill_manual(values = qualitative_colors[c(1,2)]))

# fuzzy - quality - BCI
fuzzy_quality_BCI <- FuzzyQ_plants %>% 
  mutate(evaluate_Sil = case_when(sil_width > 0.71 ~ "Strong cluster",
                                  sil_width > 0.51 ~ "Reasonable cluster",
                                  sil_width > 0.26 ~ "Weak cluster",
                                  sil_width <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "FuzzyQ quality - BCI dataset")

(fuzzy_plants_plot <- 
    FuzzyQ_plants %>% 
    left_join(BCI_tidy, by = c("Species" = "Taxa_id")) %>% 
    group_by(Sample) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
    ungroup() %>% 
    mutate(Classification = ifelse(cluster == 0, "Rare", "Common")) %>% 
    ggplot(aes(reorder(Species, -RelativeAbundance),
               RelativeAbundance, 
               col = Classification)) + 
    stat_summary() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          legend.position = "top") + 
    labs(x = "ranked plant species",
         y = "mean \U00B1 sd relative abundance (%)",
         col = "classification: ",
         title = "FuzzyQ classification of plants dataset") +
    scale_color_manual(values = qualitative_colors[c(7,3)]))

# remove this step #
fuzzy_vs_ulrb <-
  all_years_fuzzy$spp[,c(1,3)] %>% 
  as.data.frame() %>%
  mutate(Sequence = rownames(.)) %>% 
  left_join(all_years_ulrb)

# also remove this figure
grid.arrange(
  all_years_ulrb %>%
    group_by(Sample, year) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>%
    ungroup() %>% 
    ggplot(aes(reorder(Sequence, -RelativeAbundance), 
               RelativeAbundance, col = Classification)) + 
    stat_summary(alpha = 0.7) + 
    facet_grid(~year, scales = "free_x") + 
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(), 
          strip.background = element_blank(),
          legend.position = "top") + 
    labs(x = "ranked ASVs",
         y = "relative abundance \n(mean \U00B1 sd, Log10 scale)",
         col = "ulrb classification",
         tag = "a") + 
    scale_y_log10() + 
    geom_hline(yintercept = c(0.01, 0.1, 1), linetype = "dashed")
  ,
  fuzzy_vs_ulrb %>% 
    mutate(fuzzy = ifelse(cluster == 0, "Common", "Rare")) %>% 
    group_by(Sample, year) %>% 
    mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>%
    ungroup() %>% 
    ggplot(aes(reorder(Sequence, -RelativeAbundance), 
               RelativeAbundance, col = fuzzy)) + 
    stat_summary(alpha = 0.7) + 
    facet_grid(~year, scales = "free_x") + 
    scale_color_manual(values = qualitative_colors[c(7,3)]) + 
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank(), 
          strip.background = element_blank(),
          legend.position = "top") + 
    labs(x = "ranked ASVs",
         y = "relative abundance \n(mean \U00B1 sd, Log10 scale)",
         col = "FuzzyQ classification",
         tag = "b") + 
    scale_y_log10() + 
    geom_hline(yintercept = c(0.01, 0.1, 1), linetype = "dashed"),
  nrow = 2
)

### ulrb with FuzzyQ data
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
# 
gridExtra::grid.arrange(
ants_ulrb %>%
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>%
  ungroup() %>% 
  mutate(Taxa_id = factor(Taxa_id, levels = ranked_ants)) %>% 
  ggplot(aes(Taxa_id, y = Abundance, col = Classification)) +
  #geom_point(size = 1.5) +
  stat_summary() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top") + 
  labs(x = "ranked ant species",
       y = "relative abundance (%) (mean \U00B1 sd)",
       col = "classification: ",
       title = "ulrb - Ants dataset") +
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
        legend.position = "top") + 
  geom_hline(yintercept = 0)  + 
  geom_hline(yintercept = 0.5, lty = "dashed") + 
  labs(x = "ranked ant species",
       y = "Silhouette score (mean \U00B1 sd)",
       col = "Classification ",
       title = "ulrb - Ants dataset") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)]))

#
ulrb_quality_ants <- 
  ants_ulrb %>% 
  mutate(evaluate_Sil = case_when(Silhouette_scores > 0.71 ~ "Strong cluster",
                                  Silhouette_scores > 0.51 ~ "Reasonable cluster",
                                  Silhouette_scores > 0.26 ~ "Weak cluster",
                                  Silhouette_scores <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "ulrb quality - Ants dataset")


# ulrb - mosj2016-2020 (full data, similar to FuzzyQ)
gridExtra::grid.arrange(
all_years_ulrb %>%
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ggplot(aes(reorder(Sequence, -Abundance),
             RelativeAbundance,
             col = Classification)) + 
  stat_summary() + 
  scale_y_log10() + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top") + 
  labs(x = "ranked ASVs",
       y = "Relative abundance (%) (mean \U00B1 sd)",
       col = "Classification ",
       title = "ulrb - MOSJ2016-2020 dataset") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)]),
all_years_ulrb %>% 
  ggplot(aes(reorder(Sequence, -Silhouette_scores),
             Silhouette_scores,
             col = Classification)) + 
  stat_summary() +   
  theme_classic() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top") + 
  labs(x = "ranked ASVs",
       y = "Silhouette scores (mean \U00B1 sd)",
       col = "Classification ") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") + 
  facet_grid(~Classification))
  
# Count
ulrb_quality_mosj <-
all_years_ulrb %>% 
  mutate(evaluate_Sil = case_when(Silhouette_scores > 0.71 ~ "Strong cluster",
                                  Silhouette_scores > 0.51 ~ "Reasonable cluster",
                                  Silhouette_scores > 0.26 ~ "Weak cluster",
                                  Silhouette_scores <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "ulrb quality - MOSJ2016-2020 dataset")
  
  
# remove this part
ants_fuzzy_q <- antsA %>% 
  fuzzyq(sorting = TRUE)

ants_fuzzy_q_tidy <- ants_fuzzy_q$spp[, c(1,3)] %>% 
  as.data.frame() %>% 
  rename(Classification = cluster,
         Probability = Common.I) %>% 
  mutate(Species = rownames(ants_fuzzy_q$spp)) %>% 
  mutate(Classification = ifelse(Classification == 0, "Rare", "Common")) %>% 
  full_join(antsA_tidy, by = c("Species" = "Taxa_id"))

# merge fuzzyQ ants results with samples
fuzzyq_ants_plot <- 
  ants_fuzzy_q_tidy %>% 
  ggplot(aes(reorder(Species, -Abundance), Abundance, col = Classification)) + 
  stat_summary() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top") + 
  labs(x = "ranked ant species",
       y = "mean \U00B1 sd relative abundance (%)",
       col = "classification: ",
       title = "FuzzyQ classification of ants dataset") +
  scale_color_manual(values = qualitative_colors[c(7,3)])

### FuzzyQ vs ulrb

# transform to long format
BCI_tidy <- prepare_tidy_data(BCI, samples_in = "rows", sample_names = rownames(BCI))
# 
BCI_ulrb <- BCI_tidy %>% 
  define_rb()

# BCI - ulrb
gridExtra::grid.arrange(
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
                          legend.position = "top") + 
  labs(x = "ranked ant species",
       y = "relative abundance (%) (mean \U00B1 sd )",
       col = "Classification ",
       title = "ulrb - BCI dataset") +
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
        legend.position = "top") + 
  labs(x = "ranked tree species",
       y = "Silhouette score (mean \U00B1 sd)",
       col = "Classification") +
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = c(0.5), lty = "dashed"))
  
# ulrb quality - BCI
ulrb_quality_BCI <-
  BCI_ulrb %>% 
  mutate(evaluate_Sil = case_when(Silhouette_scores > 0.71 ~ "Strong cluster",
                                  Silhouette_scores > 0.51 ~ "Reasonable cluster",
                                  Silhouette_scores > 0.26 ~ "Weak cluster",
                                  Silhouette_scores <= 0.26 ~ "Potentially artificial")) %>% 
  group_by(Classification) %>% 
  count(evaluate_Sil) %>% 
  mutate(Prop = n*100/sum(n)) %>% 
  mutate(evaluate_Sil = factor(evaluate_Sil, levels = c("Strong cluster",
                                                        "Reasonable cluster",
                                                        "Weak cluster",
                                                        "Potentially artificial"))) %>% 
  ggplot(aes(Classification, Prop, fill = evaluate_Sil)) + 
  geom_col() + 
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#F0E442", "#D55E00")) + 
  theme_bw() + 
  theme(legend.position = "top",
        y = "Proportion (%)",
        panel.grid = element_blank()) + 
  geom_hline(yintercept = 50, lty = "dashed") + 
  labs(fill = "Evaluation",
       title = "ulrb quality - BCI dataset")
  

## compare clustering quality
gridExtra::grid.arrange(
ulrb_quality_BCI + 
  guides(fill = FALSE) + 
  labs(y = "Proportion (%)"),
ulrb_quality_mosj + 
  guides(fill = FALSE)+ 
  labs(y = "Proportion (%)"),
ulrb_quality_ants + 
  guides(fill = FALSE)+ 
  labs(y = "Proportion (%)"),
fuzzy_quality_BCI + 
  guides(fill = FALSE)+ 
  labs(y = "Proportion (%)"),
fuzzy_ants_quality + 
  guides(fill = FALSE)+ 
  labs(y = "Proportion (%)"),
fuzzyQ_quality_mosj + 
  guides(fill = FALSE)+ 
  labs(y = "Proportion (%)"),
ncol = 3, nrow = 2)

###


# remove this plot
plant_clustering <- 
  BCI_ulrb %>% group_by(Sample) %>% 
  mutate(Abundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(reorder(Taxa_id, -Abundance), Abundance, col = Classification)) + 
  stat_summary()+
  theme_bw()+ 
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top") + 
  labs(x = "ranked plant species",
       y = "mean \U00B1 sd relative abundance (%)",
       title = "ulrb classification of plants dataset",
       color = "classification: ") + 
  scale_color_manual(values = qualitative_colors[c(3,4,7)])

#




grid.arrange(arrangeGrob(ulrb_ants_plot, fuzzyq_ants_plot, ncol = 2),
             arrangeGrob(plant_clustering, fuzzy_plants_plot, ncol = 2)
)


### timer for FuzzyQ vs ulrb
library(microbenchmark)

## microbenchmark is commented to avoid running it by mistake (it takes a lot of RAM and time)
#fuzzyQ_vs_ulrb_bench <- microbenchmark(
#  microbiome_data_ulrb_all_AO_samples = {
#    define_rb(all_years)
#  },
#  micobiome_data_fuzzyQ_all_AO_samples = {
#    fuzzyq(all_years_ulrb_matrix, sorting = TRUE)
#  },
#  ants_data_ulrb_all_samples = {
#    define_rb(antsA_tidy)
#  },
#  ants_data_fuzzy_q_all_samples = {
#    fuzzyq(antsA[-95,], sorting = TRUE)
#  },
#  plants_data_ulrb_all_samples = {
#    define_rb(BCI_tidy)
#  },
#  plants_data_fuzzyQ_all_samples = {
#    fuzzyq(BCI, sorting = TRUE)
#  }
#)

# don't run
#load("fuzzyQvsulrbBench")

grid.arrange(
  # Microbiome data
  fuzzyQ_vs_ulrb_bench %>% 
    filter(expr %in% c("microbiome_data_ulrb_all_AO_samples",
                       "micobiome_data_fuzzyQ_all_AO_samples")) %>% 
    mutate(expr = ifelse(expr == "microbiome_data_ulrb_all_AO_samples", "ulrb", "FuzzyQ")) %>% 
    autoplot() + 
    #  scale_x_log10() + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) + 
    labs(x = "time (miliseconds)",
         y = "Approach",
         title = "Microbiome dataset",
         subtitle = "100 replications")
  ,
  # Ants data
  fuzzyQ_vs_ulrb_bench %>% 
    filter(expr %in% c("ants_data_ulrb_all_samples",
                       "ants_data_fuzzy_q_all_samples")) %>% 
    mutate(expr = ifelse(expr == "ants_data_ulrb_all_samples", "ulrb", "FuzzyQ")) %>% 
    autoplot() +
    # scale_x_log10() + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) + 
    labs(x = "time (miliseconds)",
         y = "Approach",
         title = "Ants dataset",
         subtitle = "100 replications"),
  fuzzyQ_vs_ulrb_bench %>% 
    filter(expr %in% c("plants_data_ulrb_all_samples",
                       "plants_data_fuzzyQ_all_samples")) %>% 
    mutate(expr = ifelse(expr == "plants_data_ulrb_all_samples", "ulrb", "FuzzyQ")) %>% 
    autoplot() +
    # scale_x_log10() + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) + 
    labs(x = "time (miliseconds)",
         y = "Approach",
         title = "Plants dataset (BCI)",
         subtitle = "100 replications"),
  nrow= 1
)

