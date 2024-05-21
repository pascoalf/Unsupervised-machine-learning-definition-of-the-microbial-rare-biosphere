## testing ##
short_vs_full_ulrb %>% 
  ggplot(aes(reorder(ID, -RelativeAbundance), 
             RelativeAbundance, 
             col = Classification,
             group = Sample
             )) + 
  geom_line() + 
  #geom_point(size = 1, alpha = 0.5) + 
  #stat_summary(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = c(0.01, 0.1, 1), linetype = "dashed") +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x =element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 10)) + 
  facet_grid(~Marker, scales = "free") + 
  scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
  scale_y_log10() + 
  labs(x = "ranked ASV",
       y = "relative abundance (%) \n(Log10 scale)",
       col = "classification: ")

#log_ratio = map(.x = data,~decostand(.x$Abundance, method = "rclr"))
#
all_years_ulrb %>%
  group_by(Sample, year, Classification) %>% 
  summarise(Diversity = vegan::specnumber(initialAbundance)) %>%
  ungroup() %>% 
  ggplot(aes(Classification, Diversity))+
  geom_jitter(col = "grey42", height = 0, width = 0.15) + 
  stat_summary()+ 
  stat_summary(aes(y = Diversity, group = 1), 
               fun = mean, geom = "line",
               lwd = 1, lty = "dashed") + 
  theme_bw()

