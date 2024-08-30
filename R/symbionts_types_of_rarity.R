# types of rarity in octocoral microbiome

# load source data (courtesy if Keller-Costa et al)
symbiose_tidy <- readRDS("./source_data/octocoral_microbiome.rds")


# baseline
baseline_samples <- symbiose_tidy %>% ungroup() %>% select(cleanSample) %>% distinct()


# make function to survey specific OTUs
survey_otu <- function(data, otu){
  # taxonomy
  taxonomy <- data %>% 
    filter(OTU == otu) %>% 
    select(Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
    distinct() %>% 
    mutate(taxonomy = paste(Phylum, Class, Order, Family, Genus,sep = ", ")) %>% 
    pull(taxonomy)
  #
  data %>% 
    filter(OTU == otu) %>% 
    right_join(baseline_samples) %>% 
    mutate(Classification = as.character(Classification),
           Classification = ifelse(is.na(Classification), "Absent", Classification),
           Classification = factor(Classification, 
                                   levels = c("Absent", "Rare", "Undetermined", "Abundant"))) %>%
    mutate(group = case_when(str_detect(cleanSample, "EG16") ~ "EG16",
                             str_detect(cleanSample, "EG18") ~ "EG18",
                             str_detect(cleanSample, "EG15") ~ "EG15",
                             str_detect(cleanSample, "SD") ~ "Sediment",
                             str_detect(cleanSample, "SW") ~ "Seawater",
                             str_detect(cleanSample, "LS") ~ "LS",
                             str_detect(cleanSample, "EV") ~ "EV")) %>% 
    ggplot(aes(x = cleanSample, y = Classification)) + 
    geom_point(size = 3) + 
    geom_line(aes(group = group)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 12, angle = 90),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14)) + 
    geom_vline(xintercept = c(2.5, 4.5, 6.5, 10.5, 13.5, 16.5),
               lty = "dashed") + 
    labs(x = "Sample",
         title = otu,
         subtitle = taxonomy) +
    ylim(c("Absent", "Rare", "Undetermined", "Abundant"))
}

#
grid.arrange(
survey_otu(data = symbiose_tidy,
           otu = "OTU_561"),
survey_otu(data = symbiose_tidy,
           otu = "OTU_559"),
survey_otu(data = symbiose_tidy,
           otu = "OTU_866"), 
ncol = 1)



