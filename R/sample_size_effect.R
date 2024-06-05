## effect of sample size

##
sample_size_groups <- seq(from = 6, to = 114, by = 2)

## divide by classification
ulrb_sample_size_sil.2 <- function(size){
  available_samples <- unique(all_years_rarefaction$Sample)
  set.seed(123)
  suppressWarnings(all_years_rarefaction %>%
                     filter(Sample %in% sample(available_samples, size, replace  = FALSE)) %>%
                     define_rb() %>%
                     select(Sample, Classification, Silhouette_scores) %>%
                     group_by(Sample, Classification) %>% 
                     summarise(avgSil = mean(Silhouette_scores)) %>% 
                     mutate(size = size)
                   )
}

#
ulrb_sample_size_effect.2 <- lapply(sample_size_groups, function(x) ulrb_sample_size_sil.2(size = x))

#
ulrb_sample_size_effect.2_df <- bind_rows(ulrb_sample_size_effect.2) %>% 
  mutate(method = "ulrb")

#
ulrb_sample_size_effect.2_df %>% 
  ggplot(aes(x = size, avgSil, col = Classification))+
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") + 
  scale_x_continuous(breaks = seq(from = 0, to = 110, by = 10)) + 
  theme_classic() + 
  theme(legend.position = "top") +
  labs(x = "n", 
       y = "mean (\U00B1 sd) of average Silhouette score",
       title = "ulrb performance as a function of number of samples (n)",
       subtitle = "Total reads per sample = 10000 reads") + 
  geom_vline(xintercept = 30, lty = "dashed", col = "grey41") + 
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)])+
  scale_fill_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  ylim(0, 1) + 
  geom_hline(yintercept = c(0.7, 0.5, 0.25), lty = "dashed", col = "grey41")+ 
  geom_text(data = evaluation_sil, 
            aes(y = score+0.025, x = 110, label = evaluation), col = "black") 



