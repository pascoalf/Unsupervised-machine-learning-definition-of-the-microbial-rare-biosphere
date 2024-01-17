## Ecological analysis MOSJ 2016-2020 dataset
all_years_raw <- readRDS("./source_data/AO_2016_2020_silva.rds")

## Clean taxonomy
all_years <- all_years_raw %>% 
  filter(!is.na(Kingdom),
         !is.na(Phylum)) %>%
  mutate(Phylum = ifelse(Phylum == "Cyanobacteria/Chloroplast", 
                         "Cyanobacteria", Phylum)) %>% 
  filter(Order != "Chloroplast")

## Samples that pass rarefaction
samples_to_keep <- 
  all_years %>% 
  group_by(Sample) %>% 
  summarise(Total_reads = sum(Abundance)) %>% 
  ungroup() %>% 
  filter(Total_reads >= 10000) %>% 
  pull(Sample)

set.seed(123)
# will give warning, this is not a problem
all_years_rarefaction <- 
  all_years %>% 
  filter(Sample %in% samples_to_keep) %>%  ## filter samples with more than 10 000 reads
  filter(Abundance > 1) %>% 
  group_by(Sample) %>%
  nest() %>% 
  mutate(Rarefied_reads = map(.x = data, 
                              ~as.data.frame(
                                t(
                                  rrarefy(
                                    .x$Abundance, 
                                    sample = 10000))))) %>% 
  unnest(c(data,Rarefied_reads))%>% 
  rename(Rarefied_abundance = "V1")

## Rare Biosphere ##
all_years_ulrb <- 
  all_years_rarefaction %>%  
  mutate(Abundance = Rarefied_abundance, 
         Rarefied_abundance = NULL,
         Rarefaction = "Yes") %>% 
  define_rb()
#


all_years_richness <- 
  all_years_ulrb %>% 
  group_by(Sample, year, Depth, Classification) %>% 
  summarise(Species_richness = vegan::specnumber(Abundance)) %>% 
  ungroup()


### group by depth
(alpha_diversity_plot_2016_2020 <- 
    all_years_ulrb %>%
    mutate(Zone = case_when(Depth <= 200 ~ "epipelagic (0-200m)",
                            Depth <= 1000 ~ "mesopelagic (200-1000m)",
                            Depth >1000 ~ "bathypelagic (1000-4000m)")) %>%
    filter(!is.na(Zone)) %>%  ## check this
    mutate(Zone = factor(Zone, levels = c("epipelagic (0-200m)", 
                                          "mesopelagic (200-1000m)", 
                                          "bathypelagic (1000-4000m)"))) %>% 
    group_by(Sample, year, Zone, Classification) %>% 
    summarise(Diversity = vegan::specnumber(Abundance)) %>%
    ungroup() %>% 
    ggplot(aes(factor(year), Diversity))+
    geom_half_boxplot(outlier.shape = "cross", 
                      outlier.size = 0.5, 
                      outlier.color = "red",
                      aes(fill = Classification)) +
    geom_half_point(#size = 0.75, 
      aes(col = Classification)) +
    scale_color_manual(values = qualitative_colors[c(3, 4, 7)])+
    scale_fill_manual(values = qualitative_colors[c(3, 4, 7)])+
    labs(x = "year",
         y = "number of ASVs",
         color = "classification: ",
         fill = "classification: ") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          #legend.position = "top",
          panel.grid.major.x = element_line(color = "grey"),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) + 
    facet_grid(~Zone))



rare_phyla_count <- all_years_ulrb %>% 
  filter(Classification == "Rare") %>% 
  group_by(Sample, year, Depth) %>% 
  count(Kingdom, Phylum)
#
(taxonomy_plot <- 
    rare_phyla_count %>% 
    mutate(water_group = case_when(Depth <= 200 ~ "Epipelagic",
                                   Depth <= 1000 ~ "Mesopelagic",
                                   Depth > 1000 ~ "Bathypelagic")) %>% 
    mutate(water_group = factor(water_group, 
                                levels = c("Epipelagic", 
                                           "Mesopelagic", 
                                           "Bathypelagic"))) %>% 
    ggplot(aes(reorder(Phylum, n), n, col = water_group)) +
    stat_summary() + 
    scale_y_log10() +
    coord_flip() + 
    theme_bw() + 
    labs(x = "phyla (Silva database v132)",
         y = "mean (\U00B1 sd) number of ASVs (Log10 scale)",
         col = "pelagic zone: ",
         title = "Rare biosphere") + 
    facet_grid(~year) + 
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          #legend.position = "top"
    ) + 
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 9, "BuPu")[c(3,6,9)]))

grid.arrange(alpha_diversity_plot_2016_2020,
             taxonomy_plot, 
             layout_matrix = rbind(c(1,1,1),
                                   c(1,1,1),
                                   c(2,2,2),
                                   c(2,2,2),
                                   c(2,2,2)))

## Beta diversity ##
## Data preparation for beta diversity
all_years_ulrb_helper <- 
  all_years_ulrb %>% 
  ungroup() %>%   
  select(Sample, Classification, Abundance, Sequence, Depth, year) %>% 
  mutate(Sample_unique = paste(Sample, year, Classification))

all_years_env <- 
  all_years_ulrb_helper %>% 
  ungroup() %>% 
  select(Sample_unique, Depth, year, Classification) %>% 
  distinct()

all_years_env <- 
  all_years_env %>% 
  mutate(col_year = case_when(year == 2016 ~ qualitative_colors[2],
                              year == 2017 ~ qualitative_colors[3],
                              year == 2018 ~ qualitative_colors[4],
                              year == 2019 ~ qualitative_colors[5],
                              year == 2020 ~ qualitative_colors[6])) %>% 
  mutate(col_depth = case_when(Depth <= 200 ~ "#DEEBF7",
                               Depth <= 1000 ~ "#4292C6",
                               Depth > 1000 ~ "#08306B")) %>% 
  mutate(water_group = case_when(Depth <= 200 ~ "Epipelagic",
                                 Depth <= 1000 ~ "Mesopelagic",
                                 Depth > 1000 ~ "Bathypelagic")) %>% 
  mutate(detailed_depth = case_when(Depth <= 10 ~ brewer.pal(9, "Blues")[1],
                                    Depth <= 50 ~ brewer.pal(9, "Blues")[2],
                                    Depth <= 100 ~ brewer.pal(9, "Blues")[3],
                                    Depth <= 200 ~ brewer.pal(9, "Blues")[4],
                                    Depth <= 250 ~ brewer.pal(9, "Blues")[5],
                                    Depth <= 500 ~ brewer.pal(9, "Blues")[6],
                                    Depth <= 1000 ~ brewer.pal(9, "Blues")[7],
                                    Depth <= 2000 ~ brewer.pal(9, "Blues")[8],
                                    Depth > 2000 ~ brewer.pal(9, "Blues")[9]))

#
all_years_ulrb_helper <- all_years_ulrb_helper %>% select(-Classification,
                                                          -Depth,
                                                          -year, 
                                                          -Sample)
all_years_wide <- 
  all_years_ulrb_helper %>%  
  tidyr::pivot_wider(names_from = "Sequence",
                     values_from = "Abundance",
                     values_fn = list(count = list)) %>% 
  unchop(everything())

# turn NAs into zeros
all_years_wide[is.na(all_years_wide)] <- 0

#
rare_biosphere <- all_years_wide %>% 
  filter(str_detect(Sample_unique, "Rare"))

rare_nMDS <- metaMDS(rare_biosphere[, -1])

rare_env <- all_years_env %>% 
  filter(Sample_unique %in% rare_biosphere$Sample_unique)

plot(rare_nMDS, 
     display = "sites", 
     type = "p", 
     main = "Rare biosphere")

points(rare_nMDS,
       bg = rare_env$col_depth,
       pch = 21, 
       col = "grey", 
       cex = 2)

with(rare_env,
     ordiellipse(rare_nMDS,
                 year,
                 #draw = "polygon",
                 lty = "dashed",
                 col = "black", 
                 label = TRUE))
with(rare_env,
     ordispider(rare_nMDS,
                year,
                draw = "polygon",
                lty = "dashed",
                col = "lightgrey", 
                label = TRUE))

## Permanova to check years
permanova_results <- adonis2(rare_biosphere[, -1] ~ year*water_group, data = rare_env)
#store in csv file
#write.csv(broom::tidy(permanova_results), "./permanova_results.csv")

# betadisper to check
set.seed(123)
permutest(
  betadisper(
    vegdist(rare_biosphere[, -1]),
    group = rare_env$year,
    type = "centroid"),
  permutations = 999)

#
permutest(
  betadisper(
    vegdist(rare_biosphere[, -1]),
    group = rare_env$water_group,
    type = "centroid"),
  permutations = 999)

par(mfrow = c(1,2))
plot(betadisper(
  vegdist(rare_biosphere[, -1]),
  group = rare_env$year,
  type = "centroid"),
  main = "Verify homogeneity of variance per year")
#
plot(betadisper(
  vegdist(rare_biosphere[, -1]),
  group = rare_env$water_group,
  type = "centroid"),
  main = "Verify homogeneity of variance per pelagic layer")

## Top phyla RAC

taxonomy_count <- 
  all_years_ulrb %>% 
  group_by(Sample, year, Depth) %>% 
  count(Kingdom, Phylum, Class)
#


phyla_count <- all_years_ulrb %>% 
  group_by(Sample, year, Depth) %>% 
  count(Kingdom, Phylum)

##
proteobacteria_rac <- 
  all_years_ulrb %>% 
  filter(Phylum == "Proteobacteria") %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ggplot(aes(reorder(Sequence, -RelativeAbundance), RelativeAbundance, col = Classification)) + 
  stat_summary() +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
  scale_y_log10() + 
  labs(title = "Rank Abundance Curve for Proteobacteria",
       x = "Ranked ASVs",
       y = "Relative abundance % (Log10 scale)") +
  geom_hline(yintercept = c(0.01, 0.1, 1, 10), 
             col = "grey", lty = "dashed")
#

all_rac <- 
  all_years_ulrb %>% 
  group_by(Sample) %>% 
  mutate(RelativeAbundance = Abundance*100/sum(Abundance)) %>% 
  ggplot(aes(reorder(Sequence, -RelativeAbundance), RelativeAbundance, col = Classification)) + 
  stat_summary() +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
  scale_y_log10() + 
  labs(x = "ranked ASVs",
       y = "relative abundance (%) \n(mean \U00B1 sd - Log10 scale)",
       title = "All phyla",
       tag = "a") +
  geom_hline(yintercept = c(0.01, 0.1, 1, 10), col = "grey", lty = "dashed")
#

## RAC for other phyla ##
specific_rac <- function(x){
  # relative abundance of members of specified phylum
  all_years_ulrb %>% 
    group_by(Sample) %>% 
    mutate(RelativeAbundance = round(Abundance*100/sum(Abundance), 2)) %>% 
    filter(Phylum == x) %>%
    ggplot(aes(reorder(Sequence, -RelativeAbundance), RelativeAbundance, col = Classification)) + 
    stat_summary() +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "top",
          axis.title.y = element_text(size = 10)) +
    scale_color_manual(values = qualitative_colors[c(3,4,7)]) + 
    scale_y_log10() + 
    labs(title = x,
         x = "ranked ASVs",
         y = "relative abundance (%) \n(mean \U00B1 sd - Log10 scale)") +
    geom_hline(yintercept = c(0.01, 0.1, 1, 10), col = "grey", lty = "dashed")  
}

# phyla with more ASVS overall
all_years_ulrb %>% 
  group_by(Phylum) %>% 
  summarise(Total = n()) %>% 
  arrange(desc(Total))

proteobacteria_rac <- specific_rac("Proteobacteria") + labs(tag = "b")
bacteroidota_rac <- specific_rac("Bacteroidota") + labs(tag = "c")
#cyanobacteria_rac <- specific_rac("Cyanobacteria")
planctomycetota_rac <- specific_rac("Planctomycetota")+ labs(tag = "e")
verrucomicrobiota_rac <- specific_rac("Verrucomicrobiota")+ labs(tag = "d")
#chloroflexi_rac <- specific_rac("Chloroflexi")+ labs(tag = "b")
#

grid.arrange(
  arrangeGrob(all_rac, ncol =1, nrow = 1),
  arrangeGrob(
    proteobacteria_rac,
    bacteroidota_rac,
    verrucomicrobiota_rac,
    planctomycetota_rac,
    ncol = 2, nrow = 2)
)





