## FuzzyQ vs ulrb ##
data(antsA)

FQAnts <- fuzzyq(antsA[-95,], sorting = TRUE)

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


fuzzy_vs_ulrb <-
  all_years_fuzzy$spp[,c(1,3)] %>% 
  as.data.frame() %>%
  mutate(Sequence = rownames(.)) %>% 
  left_join(all_years_ulrb)


fuzzy_vs_ulrb %>% 
  ggplot(aes(reorder(Sequence, -Abundance), Abundance))+
  stat_summary(aes(col = as.factor(cluster))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        legend.position = "top")+
  scale_y_log10()


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


ants_ulrb <- define_rb(antsA_tidy)

ulrb_ants_plot <- 
  ants_ulrb %>% 
  ggplot(aes(reorder(Taxa_id, -Abundance), y = Abundance, col = Classification)) +
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
       title = "ulrb classification of ants dataset") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)])

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
# load BCI dataset
data(BCI)
# transform to long format
BCI_tidy <- prepare_tidy_data(BCI, samples_in = "rows", sample_names = rownames(BCI))
# 
BCI_ulrb <- BCI_tidy %>% 
  define_rb()
#
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
FuzzyQ_plants <- fuzzyq(BCI, sorting = TRUE)
FuzzyQ_plants <-  FuzzyQ_plants$spp[,c(1,3)] %>% 
  as.data.frame() %>% mutate(Species = rownames(.))

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

