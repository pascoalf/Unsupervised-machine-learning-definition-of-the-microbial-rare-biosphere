## RACs with lines 
nice_ulrb_rarefied %>% 
  group_by(Type, Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>%
  mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
  group_by(Type, Group) %>% 
  arrange(desc(Abundance)) %>% 
  mutate(uniqueRank = row_number()) %>%
  ungroup() %>% 
  ggplot(aes(x = uniqueRank, 
             y = RelativeAbundance, 
             col = Classification)) +
  geom_point(alpha = 0.75, size = 2)+
  geom_line(aes(group = Group), col = "grey72") +
  facet_grid(~Type, scales = "free") + 
  scale_color_manual(values = qualitative_colors[c(3,5,7)]) + 
  theme_bw() + 
  scale_y_log10() + 
  geom_hline(yintercept = c(1, 0.1, 0.01), lty = "dashed") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12)) + 
  labs(x = "Ranked phylogenetic units",
       y = "Relative abundance (%) \n(Log10 scale)",
       col = "Classification: ",
       tag = "a")



# silhouette scores
nice_ulrb_rarefied %>% 
  mutate(Type = factor(Type, levels = c("ASV", "OTU", "mOTU"))) %>% 
  mutate(Group = paste(Sample, Classification, sep = "_")) %>% 
  group_by(Type, Group) %>% 
  arrange(desc(Silhouette_scores)) %>% 
  mutate(uniqueRank = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(uniqueRank, 
             Silhouette_scores, col = Classification)) + 
  geom_line(aes(group = Group), col = "grey82") +
  geom_point(alpha = 0.5, size = 2) +
  facet_wrap(~Type, scales = "free_x") + 
  scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12)) +
  guides(color = "none") + 
  labs(y = "Silhouette score \n(Log10 scale)",
       x = "Ranked phylogenetic units",
       tag = "b") + 
  geom_hline(yintercept = c(0)) + 
  geom_hline(yintercept = c(-0.5, 0.5), lty = "dashed", color = "grey41")


# end
