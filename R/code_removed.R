# code removed

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


# 
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