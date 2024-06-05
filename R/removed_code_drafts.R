## removed code


ulrb_sample_size_sil <- function(size){
  set.seed(123)
  suppressWarnings(all_years_rarefaction %>%
                     filter(Sample %in% sample(mosj_samples, size, replace  = FALSE)) %>%
                     define_rb() %>%
                     pull(Silhouette_scores) %>% 
                     mean())
  
}

#
sample_size_groups <- seq(from = 10, to = 100, by = 10)

ulrb_sample_size_effect <- sapply(sample_size_groups, function(x) ulrb_sample_size_sil(size = x))

ulrb_sample_size_effect_df <- ulrb_sample_size_effect %>% 
  as.data.frame() %>% 
  rename(Sil = ".") %>% 
  mutate(sample_size = sample_size_groups,
         method = "ulrb")

#
ulrb_sample_size_effect_df %>% 
  ggplot(aes(sample_size, Sil)) + 
  geom_point() + 
  geom_line(aes(group = 1)) + 
  ylim(c(0,1))


## FuzzyQ 

fuzzyq_sample_size <- function(size){
  # generate random group of samples
  set.seed(123); random_samples_numbers <- sample(seq(1,114), size)
  # make subsampled matrix
  sub_matrix <- all_years_ulrb_matrix[random_samples_numbers, ]
  # remove NAs
  sub_matrix[is.na(sub_matrix)] <- 0
  # apply FuzzyQ
  set.seed(123); all_years_fuzzy <- fuzzyq(sub_matrix, sorting = TRUE)
  # return average silhouette score
  all_years_fuzzy$global[3]
}

# apply to all sample sizes

fuzzyq_sample_size_effect <- sapply(sample_size_groups, function(x) fuzzyq_sample_size(size = x))

fuzzyq_sample_size_effect_df <- fuzzyq_sample_size_effect %>% 
  as.data.frame() %>% 
  rename(Sil = ".") %>% 
  mutate(sample_size = sample_size_groups,
         method = "FuzzyQ")

#
fuzzyq_sample_size_effect_df %>% 
  ggplot(aes(sample_size, Sil)) + 
  geom_point()+
  geom_line(aes(group = 1))

##
fuzzyq_sample_size_effect_df %>% 
  rbind(ulrb_sample_size_effect_df) %>% 
  ggplot(aes(sample_size, Sil, col = method)) + 
  geom_point()+
  geom_line(aes(group = method))+
  theme_classic() + 
  theme(legend.position = "top") + 
  labs(y = "average Silhouette score",
       x = "n") + 
  scale_x_continuous(breaks = sample_size_groups)

#

## same for FuzzyQ

fuzzyq_sample_size.2 <- function(size){
  # generate random group of samples
  set.seed(123); random_samples_numbers <- sample(seq(1,114), size)
  # make subsampled matrix
  sub_matrix <- all_years_ulrb_matrix[random_samples_numbers, ]
  # remove NAs
  sub_matrix[is.na(sub_matrix)] <- 0
  # apply FuzzyQ
  set.seed(123); all_years_fuzzy <- fuzzyq(sub_matrix, sorting = TRUE)
  # return average silhouette score per classification
  all_years_fuzzy$spp[,c(1:3)] %>% 
    mutate(Classification = ifelse(cluster == 0, "Rare", "Common")) %>% 
    group_by(Classification) %>% 
    summarise(avgSil = mean(sil_width)) %>% 
    mutate(size = size)
}


# apply to all combinations
fuzzyq_sample_size_effect.2 <- lapply(sample_size_groups, function(x) fuzzyq_sample_size.2(size = x))

#
fuzzyq_sample_size_effect.2_df <- bind_rows(fuzzyq_sample_size_effect.2) %>% 
  mutate(method = "FuzzyQ",
         Sample = NA)

#
fuzzyq_sample_size_effect.2_df %>% 
  ggplot(aes(x = size, avgSil, col = Classification))+
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line",
               lwd = 1, lty = "dashed") + 
  scale_x_continuous(breaks = sample_size_groups) + 
  theme_classic() + 
  theme(legend.position = "top") +
  labs(x = "n", 
       y = "mean (\U00B1 sd) of average Silhouette score") + 
  geom_vline(xintercept = 30, lty = "dashed", col = "grey41") + 
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)])+
  scale_fill_manual(values = qualitative_colors[c(3, 4, 7)])


## combine them
fuzzyq_sample_size_effect.2_df %>% 
  rbind(ulrb_sample_size_effect.2_df) %>% 
  mutate(Classification = factor(Classification, levels = c("Rare", "Undetermined", "Abundant", "Common"))) %>% 
  ggplot(aes(x = size, avgSil, col = Classification))+
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line",
               lwd = 1, lty = "dashed") + 
  scale_x_continuous(breaks = sample_size_groups) + 
  theme_classic() + 
  theme(legend.position = "top") +
  labs(x = "n", 
       y = "mean (\U00B1 sd) of average Silhouette score") + 
  geom_vline(xintercept = 30) + 
  geom_hline(yintercept = c(0.7, 0.5, 0.25), lty = "dashed", col = "grey41")+ 
  scale_color_manual(values = qualitative_colors[c(3, 4, 6, 7)])+
  scale_fill_manual(values = qualitative_colors[c(3, 4, 6, 7)]) + 
  facet_grid(~method) 
