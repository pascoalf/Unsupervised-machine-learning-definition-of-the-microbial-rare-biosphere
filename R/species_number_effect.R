## effect of number of species

#
set.seed(123); test_1000_species <- all_years_rarefaction %>% 
  group_by(Sample) %>% 
  slice_sample(n = 1000) %>% 
  define_rb() %>% 
  group_by(Sample, Classification) %>% 
  summarise(avgSil = mean(Silhouette_scores))

#
test_species_size <- function(size){
  all_years_rarefaction %>% 
    group_by(Sample) %>% 
    slice_sample(n = size) %>% 
    define_rb() %>% 
    group_by(Sample, Classification) %>% 
    summarise(avgSil = mean(Silhouette_scores)) %>% 
    mutate(size = size)
}

# species size to test
species_size_groups <- seq(from = 100, to = 4000, by = 100)

# calculate ulrb for all combinations
all_species_test <- lapply(species_size_groups, function(size) test_species_size(size = size))


all_species_test_df <- bind_rows(all_species_test)

# make labels for silhoutte scores
evaluation_sil <- data.frame(score = c(0.26, 0.51, 0.71),
                             evaluation = c("Weak cluster",
                                            "Reasonable cluster",
                                            "Strong cluster"),
                             Classification = "Rare") # placeholder for label

#
all_species_test_df %>% 
  ggplot(aes(size, avgSil, col = Classification)) + 
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") +
  ylim(0,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  theme_classic() + 
  labs(x = "number of species",
       y = "mean (\U00B1 sd) of average Silhouette score",
       title = "ulrb performance as a function of number of species",
       subtitle = "total reads per sample = 50000 reads\nn = 34 samples")+
  theme(legend.position = "top") + 
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  geom_text(data = evaluation_sil, aes(y = score+0.01, x = 500, label = evaluation), col = "black") + 
  scale_x_continuous(breaks = seq(from = 100, to = 4000, by = 300))
  
  

