# ulrb vs dataset size 

grid.arrange(
ulrb_vs_samples_plot + 
  labs(subtitle = NULL, tag = "a") + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)),
ulrb_vs_species_plot + 
  labs(subtitle = NULL,
       tag = "b") + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  guides(col = FALSE),
ulrb_vs_seq_power_plot + 
  labs(subtitle = NULL,
       tag = "c") + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  guides(col = FALSE),
layout_matrix = rbind(c(1, 1),
                      c(2, 3)))


## fuzzyQ
grid.arrange(
fuzzyq_vs_sample_size + 
  labs(subtitle = NULL, tag = "a",
       col = "Classification: ") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)),
fuzzy_vs_species_number + 
  labs(subtitle = NULL, tag = "b") + 
  guides(col = FALSE) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)),
fuzzy_vs_seq_power + 
  labs(subtitle = NULL,
       tag = "c") + 
  guides(col = FALSE) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)),
layout_matrix = rbind(c(1, 1),
                      c(2, 3)))
