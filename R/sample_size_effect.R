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
ulrb_vs_samples_plot <- ulrb_sample_size_effect.2_df %>% 
  ggplot(aes(x = size, avgSil, col = Classification))+
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") + 
  scale_x_continuous(breaks = seq(from = 0, to = 110, by = 10)) + 
  theme_classic() + 
  theme(legend.position = "top") +
  labs(x = "Number of samples (n)", 
       y = "Mean (\U00B1 sd) of \naverage Silhouette score",
       col = "Classification: ",
       #title = "ulrb performance as a function of number of samples (n)",
       subtitle = "Total reads per sample = 10000 reads") + 
  geom_vline(xintercept = 30, lty = "dashed", col = "grey41") + 
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)])+
  scale_fill_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  ylim(0, 1) + 
  geom_hline(yintercept = c(0.7, 0.5, 0.25), lty = "dashed", col = "grey41")+ 
  geom_text(data = evaluation_sil, ## data frane nade in "species_number_effect.R"
            aes(y = score+0.025, x = 110, label = evaluation), col = "black") 


## Same analysis for FuzzyQ

# change col names for easier data handling
colnames(all_years_ulrb_matrix) <- paste("ASV", seq_along(1:length(colnames(all_years_ulrb_matrix))))
#


fuzzyq_sample_size <- function(x, size){
  # subsample n samples
  set.seed(123); subsampled_data <- x[sample(1:114, size = size, replace = FALSE),]
  # replace NAs with zeros
  subsampled_data[is.na(subsampled_data)] <- 0
  #
  set.seed(123); fuzzyq_cluster <- fuzzyq(subsampled_data, sorting = TRUE)
  # return data frame with relevant metrics
  fuzzyq_cluster$spp %>% 
    as.data.frame() %>% 
    mutate(Species = rownames(.),
           Classification = ifelse(cluster == 0, "Rare", "Common")) %>% 
    group_by(Classification) %>% 
    summarise(avgSil = mean(sil_width),
              size = size)
}
#

# apply to all groups
all_sample_groups_fuzzyq <- lapply(
  sample_size_groups, 
  function(size) fuzzyq_sample_size(x = all_years_ulrb_matrix, size = size))

# silhouette labels
evaluation_sil_extra <- evaluation_sil
evaluation_sil_extra[4,] <- c(0, "Potentially artificial", "Rare")
evaluation_sil_extra$score <- as.numeric(evaluation_sil_extra$score)
#
(fuzzyq_vs_sample_size <- 
  bind_rows(all_sample_groups_fuzzyq) %>% 
  mutate(Classification = factor(Classification, levels = c("Rare", "Common"))) %>% 
  ggplot(aes(x = size, y = avgSil, col = Classification)) + 
  geom_point() +
  geom_line(aes(group = Classification)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 110, by = 10)) + 
  theme_classic() + 
  theme(legend.position = "top") +
  labs(x = "n", 
       y = "average Silhouette score",
  #     title = "fuzzyQ performance as a function of number of samples (n)",
       subtitle = "Total reads per sample = 10000 reads") + 
  geom_vline(xintercept = 30, lty = "dashed", col = "grey41") + 
  scale_color_manual(values = qualitative_colors[c(3, 7)]) +
  ylim(-0.5, 1) + 
  geom_hline(yintercept = c(0.7, 0.5, 0.25), lty = "dashed", col = "grey41")+ 
  geom_text(data = evaluation_sil_extra, ## data frane nade in "species_number_effect.R"
            aes(y = score+0.025, x = 110, label = evaluation), col = "black") )

