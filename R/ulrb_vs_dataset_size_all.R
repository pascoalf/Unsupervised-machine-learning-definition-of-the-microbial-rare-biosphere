# ulrb vs dataset size 

grid.arrange(
ulrb_vs_samples_plot + 
  labs(subtitle = NULL, tag = "a"),
ulrb_vs_species_plot + 
  labs(subtitle = NULL,
       tag = "b") + 
  guides(col = FALSE),
ulrb_vs_seq_power_plot + 
  labs(subtitle = NULL,
       tag = "c")+ 
  guides(col = FALSE),
layout_matrix = rbind(c(1, 1),
                      c(2, 3)))


## fuzzyQ
grid.arrange(
fuzzyq_vs_sample_size + 
  labs(subtitle = NULL, tag = "a"),
fuzzy_vs_species_number + 
  labs(subtitle = NULL, tag = "b") + 
  guides(col = FALSE),
fuzzy_vs_seq_power + 
  labs(subtitle = NULL,
       tag = "c") + 
  guides(col = FALSE),
layout_matrix = rbind(c(1, 1),
                      c(2, 3)))