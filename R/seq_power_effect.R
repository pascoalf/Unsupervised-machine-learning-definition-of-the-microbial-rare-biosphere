## seq power effect

nested_rarefaction <- function(size){
  #
  keep_samples <- 
    all_years %>% 
    group_by(Sample) %>% 
    summarise(Total_reads = sum(Abundance)) %>% 
    ungroup() %>% 
    filter(Total_reads >= size) %>% 
    pull(Sample)
  
  # will give warning, this is not a problem
  set.seed(123); suppressWarnings(
    all_years %>% 
    filter(Sample %in% keep_samples) %>%
    filter(Abundance > 1) %>% 
    group_by(Sample) %>%
    nest() %>% 
    mutate(Rarefied_reads = map(.x = data, 
                                ~as.data.frame(
                                  t(
                                    rrarefy(
                                      .x$Abundance, 
                                      sample = size))))) %>% 
    unnest(c(data,Rarefied_reads)) %>% 
    rename(Rarefied_abundance = "V1")  %>% 
    define_rb() %>% 
    group_by(Sample, Classification) %>% 
    summarise(avgSil = mean(Silhouette_scores)) %>% 
    mutate(size = size,
           n = length(keep_samples))
    )  
}


# total reads to test

reads_groups <- seq(from = 1000, to = 100000, by = 3000)

ulrb_vs_reads <- lapply(reads_groups, function(size) nested_rarefaction(size = size))


ulrb_vs_reads_df <- bind_rows(ulrb_vs_reads)


gridExtra::grid.arrange(
ulrb_vs_reads_df %>% 
  ggplot(aes(size, avgSil, col = Classification)) + 
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") + 
  theme_classic() + 
  theme(legend.position = "top") + 
  ylim(0,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  geom_text(data = evaluation_sil, # this object was made in "R/species_number_effect.R"
            aes(y = score+0.025, x = 10000, label = evaluation), col = "black"),
ulrb_vs_reads_df %>% 
  ggplot(aes(size, n)) + 
  geom_point() + 
  geom_hline(yintercept = 30, lty = "dashed") + 
  theme_classic() + 
  labs(x = "number of reads per sample",
       y = "n",
        title = "Reads per sample vs number of samples"), 
ncol = 2
)

# 2 variables change at the same time!
gridExtra::grid.arrange(
  ulrb_vs_reads_df %>% 
    ggplot(aes(size, avgSil, col = Classification)) + 
    stat_summary() + 
    stat_summary(aes(y = avgSil, group = Classification, 
                     color = Classification), 
                 fun = mean, geom = "line") + 
    theme_classic() + 
    theme(legend.position = "top") + 
    ylim(0,1) + 
    geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
    geom_text(data = evaluation_sil, # this object was made in "R/species_number_effect.R"
              aes(y = score+0.025, x = 10000, label = evaluation), col = "black") + 
    labs(y = "mean (\U00B1 sd) of average Silhouette score",
         x = "number of reads per sample"),
  ulrb_vs_reads_df %>% 
    ggplot(aes(n, avgSil, col = Classification)) + 
    stat_summary() + 
    stat_summary(aes(y = avgSil, group = Classification, 
                     color = Classification), 
                 fun = mean, geom = "line") + 
#    geom_hline(yintercept = 30, lty = "dashed") + 
    theme_classic() + 
    labs(y = "mean (\U00B1 sd) of average Silhouette score",
         x = "number of samples (n)") + 
  theme(legend.position = "top"), 
  ncol = 2
)


## try a different approach (confounding variables)


# Lock the number of samples
# use 50000 reads max, because we can keep 34 samples

hq_samples <- all_years %>% 
  group_by(Sample) %>% 
  summarise(Total_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  filter(Total_reads >= 50000) %>% 
  pull(Sample) 
  
# 
nested_rarefaction.2 <- function(size){
  if(size > 50000){
    stop("Maximum size is 50000")
  }
  #
  # will give warning, this is not a problem
  set.seed(123); suppressWarnings(
    all_years %>% 
      filter(Sample %in% hq_samples) %>%
      filter(Abundance > 1) %>% 
      group_by(Sample) %>%
      nest() %>% 
      mutate(Rarefied_reads = map(.x = data, 
                                  ~as.data.frame(
                                    t(
                                      rrarefy(
                                        .x$Abundance, 
                                        sample = size))))) %>% 
      unnest(c(data,Rarefied_reads)) %>% 
      rename(Rarefied_abundance = "V1")  %>% 
      define_rb() %>% 
      group_by(Sample, Classification) %>% 
      summarise(avgSil = mean(Silhouette_scores)) %>% 
      mutate(size = size,
             n = length(hq_samples))
  )  
}
#
reads_groups.2 <- seq(from = 1000, to = 50000, by = 1000)
#
ulrb_vs_reads.2 <- lapply(reads_groups.2, function(size) nested_rarefaction.2(size = size))
# 
ulrb_vs_reads_df.2 <- bind_rows(ulrb_vs_reads.2)


#
ulrb_vs_reads_df.2 %>% 
  ggplot(aes(size, avgSil, col = Classification)) + 
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") + 
  theme_classic() + 
  theme(legend.position = "top") + 
  ylim(0,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  geom_text(data = evaluation_sil, # this object was made in "R/species_number_effect.R"
            aes(y = score+0.025, x = 3500, label = evaluation), col = "black") + 
  labs(y = "mean (\U00B1 sd) of average Silhouette score",
       x = "number of reads per sample",
       title = "ulrb performance as a function of number of reads per sample",
       subtitle = "n = 34 samples") +
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)])


## fuzzyQ
multiple_rarefaction_fuzzyQ <- function(x, size){
  if(size > 50000){
    stop("Maximum size is 50000")
  }
set.seed(123); x_rarefied <- rrarefy(x, sample = size)

# remove missing values
x_rarefied[is.na(x_rarefied)] <- 0

# apply fuzzyQ
fuzzyq_cluster <- fuzzyq(x_rarefied, sorting = TRUE)

#
fuzzyq_cluster$spp %>% 
  as.data.frame() %>% 
  mutate(Species = rownames(.),
         Classification = ifelse(cluster == 0, "Rare", "Common")) %>% 
  group_by(Classification) %>% 
  summarise(avgSil = mean(sil_width),
            size = size)
}

# run fuzzyQ for all seq. powers
fuzzyQ_all_seq_power <- lapply(reads_groups.2, function(size){
  multiple_rarefaction_fuzzyQ(x = mosj_matrix_selected, size = size)})

# plot performance vs seq power
fuzzyQ_all_seq_power %>%
  bind_rows() %>% 
  ggplot(aes(size, avgSil, col = Classification)) + 
  stat_summary() + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") + 
  theme_classic() + 
  theme(legend.position = "top") + 
  ylim(0,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  geom_text(data = evaluation_sil, # this object was made in "R/species_number_effect.R"
            aes(y = score+0.025, x = 3500, label = evaluation), col = "black") + 
  labs(y = "mean (\U00B1 sd) of average Silhouette score",
       x = "number of reads per sample",
       title = "fuzzyQ performance as a function of number of reads per sample",
       subtitle = "n = 34 samples") +
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)])
