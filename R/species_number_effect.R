## effect of number of species
# (example for a single test)
#set.seed(123); test_1000_species <- all_years_rarefaction %>% 
#  group_by(Sample) %>% 
#  slice_sample(n = 1000) %>% 
#  define_rb() %>% 
#  group_by(Sample, Classification) %>% 
#  summarise(avgSil = mean(Silhouette_scores))

#
test_species_size <- function(size){
  set.seed(123)
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
ulrb_vs_species_plot <- all_species_test_df %>% 
  ggplot(aes(size, avgSil, col = Classification, fill = Classification)) + 
  geom_point(alpha = 0.1) + 
  stat_summary(size = 1, shape = 21, col = "black") + 
  stat_summary(aes(y = avgSil, group = Classification, 
                   color = Classification), 
               fun = mean, geom = "line") +
  #ylim(0,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  theme_classic() + 
  labs(x = "Number of ASVs",
       y = "Mean (\U00B1 sd) of \naverage Silhouette score",
       #title = "ulrb performance as a function of number of species",
      # subtitle = "total reads per sample = 50000 reads\nn = 34 samples"
      )+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        panel.background = element_rect(fill = "grey90")) + 
  scale_color_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  scale_fill_manual(values = qualitative_colors[c(3, 4, 7)]) + 
  #geom_label(data = evaluation_sil, 
   #          aes(y = score+0.025, x = 3000, label = evaluation), 
    #         col = "black", fill = "white") + 
  scale_x_continuous(breaks = seq(from = 100, to = 4000, by = 300))
  

## FuzzyQ

# Prepare microbiome data for fuzzyQ
mosj_matrix <- all_years %>% 
  ungroup() %>%   
  select(Sample, Abundance, Sequence, Depth, year) %>% 
  mutate(Sample_unique = paste(Sample, year)) %>%
  select(Sample_unique, Abundance, Sequence) %>% 
  tidyr::pivot_wider(names_from = "Sequence",
                     values_from = "Abundance",
                     values_fn = list(count= list)) %>% 
  print() %>% 
  unchop(everything())


rownames(mosj_matrix) <- mosj_matrix$Sample_unique
mosj_matrix$Sample_unique <- NULL
mosj_matrix <- as.matrix(mosj_matrix)

# rarefy to 50000 reads, selected 34 samples
mosj_matrix_selected <- mosj_matrix[rowSums(mosj_matrix, na.rm = TRUE) >= 50000,]
# remove NAs
mosj_matrix_selected[is.na(mosj_matrix_selected)] <- 0
mosj_matrix_rarefied_50000 <- rrarefy(mosj_matrix_selected, sample = 50000)

#
fuzzyq_species_size <- function(x, size){
  # subsample n samples
  set.seed(123); subsampled_data <- x[, sample(1:length(colnames(x)), 
                                               size = size, replace = FALSE)]
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

# all groups of species size

fuzzyq_all_species <- lapply(species_size_groups, 
                             function(size) fuzzyq_species_size(size = size, 
                                                                x = mosj_matrix_rarefied_50000))
##
fuzzyq_all_species.df <- bind_rows(fuzzyq_all_species)

#  
fuzzy_vs_species_number <- 
fuzzyq_all_species.df %>% 
  mutate(Classification = factor(Classification, levels = c("Rare", "Common"))) %>% 
  ggplot(aes(size, avgSil, col = Classification)) + 
  geom_point() + 
  geom_line(aes(group = Classification)) + 
  ylim(-0.5,1) + 
  geom_hline(yintercept = c(0.25, 0.5, 0.7), lty = "dashed") + 
  geom_vline(xintercept = 700, lty = "dashed") + 
  theme_classic() + 
  labs(x = "number of ASVs",
       y = "average Silhouette score",
#       title = "fuzzyQ performance as a function of number of species",
       #subtitle = "total reads per sample = 50000 reads\nn = 34 samples"
) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90)) + 
  scale_color_manual(values = qualitative_colors[c(3, 7)]) + 
  geom_text(data = evaluation_sil_extra, aes(y = score+0.075, x = 3500, label = evaluation), col = "black") + 
  scale_x_continuous(breaks = seq(from = 100, to = 4000, by = 300))
 

## extras
test_10 <- all_years_rarefaction %>% 
   group_by(Sample) %>% 
  slice_sample(n = 10)

test_10 %>% 
  define_rb() %>% 
  plot_ulrb_clustering(taxa_col = "Sequence", log_scaled = TRUE, sample_id = "m18_33", plot_all = FALSE)


## two examples for supplementary data
test_n <- function(n){
  all_years_rarefaction %>% 
    group_by(Sample) %>% 
    slice_sample(n = n) %>% 
    mutate(size = paste(n, "ASVs"))
  }
#
test_n(100) %>% 
 rbind(test_n(1000)) %>% 
  rbind(test_n(3000)) %>% 
  define_rb() %>% 
plot_ulrb_clustering(taxa_col = "Sequence", log_scaled = TRUE, sample_id = "m18_33", plot_all = FALSE) + 
  facet_grid(~size) + 
  labs(x = "Ranked ASVs")

